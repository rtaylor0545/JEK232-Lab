## GSE167363 Sepsis patients scRNAseq 0h ##

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(celldex)
library(SingleR)

HC1 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/HC1/")
HC1_so <- CreateSeuratObject(counts = HC1, project = "HC1")
HC1_so[["percent.mt"]] = PercentageFeatureSet(HC1_so, pattern = "^MT.")
HC2 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/HC2/")
HC2_so <- CreateSeuratObject(counts = HC2, project = "HC2")
HC2_so[["percent.mt"]] = PercentageFeatureSet(HC2_so, pattern = "^MT.")

GD1 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/SS1/")
GD1_so <- CreateSeuratObject(counts = GD1, project = "GD1")
GD1_so[["percent.mt"]] = PercentageFeatureSet(GD1_so, pattern = "^MT.")
GD2 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/SS2/")
GD2_so <- CreateSeuratObject(counts = GD2, project = "GD2")
GD2_so[["percent.mt"]] = PercentageFeatureSet(GD2_so, pattern = "^MT.")
GD3 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/SS3/")
GD3_so <- CreateSeuratObject(counts = GD3, project = "GD3")
GD3_so[["percent.mt"]] = PercentageFeatureSet(GD3_so, pattern = "^MT.")

BD1 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/NS1/")
BD1_so <- CreateSeuratObject(counts = BD1, project = "BD1")
BD1_so[["percent.mt"]] = PercentageFeatureSet(BD1_so, pattern = "^MT.")
BD2 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/6h/NS2/")
BD2_so <- CreateSeuratObject(counts = BD2, project = "BD2")
BD2_so[["percent.mt"]] = PercentageFeatureSet(BD2_so, pattern = "^MT.")

# Merging all samples
pbmc.combined <- merge(HC1_so, y = c(HC2_so, GD1_so, GD2_so, GD3_so, BD1_so, BD2_so), 
                       add.cell.ids = c("HC1", "HC2", "GD1", "GD2", "GD3", "BD1", "BD2"), 
                       project = "GSE167363_6h")

table(pbmc.combined$orig.ident)

#Scatter 플랏
plot0 <- VlnPlot(pbmc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width = 12, repr.plot.height = 6)
plot0 
plot1 + plot2

rm(HC2_so)

pbmc.combined2 <- subset(pbmc.combined, subset = nFeature_RNA > 50 & nFeature_RNA < 5000 & percent.mt < 30)
pbmc.combined2 <- NormalizeData(pbmc.combined2, verbose = FALSE)
pbmc.combined2 <- FindVariableFeatures(pbmc.combined2, verbose = FALSE)
pbmc.combined2 <- ScaleData(pbmc.combined2, verbose = FALSE)
pbmc.combined2 <- RunPCA(pbmc.combined2, verbose = FALSE)

saveRDS(pbmc.combined, file = "D:/HP_rdata/GSE167363/6h_Sepsis_merging_seurat_data(prototype).rds")
saveRDS(pbmc.combined2, file = "D:/HP_rdata/GSE167363/6h_Sepsis_merging_seurat_data(filter).rds")
write.csv(pbmc.combined@active.ident, file = "d:/HP_rdata/GSE167363/6h_cell_population(proto).csv")
write.csv(pbmc.combined2@active.ident, file = "d:/HP_rdata/GSE167363/6h_cell_population(filter).csv")

############# Analysis ###############
pbmc.combined2 <- readRDS(file = "D:/HP_rdata/GSE167363/6h_Sepsis_merging_seurat_data(filter).rds")

plot <- ElbowPlot(pbmc.combined2)
options(repr.plot.width = 6, repr.plot.height = 6)
plot

pbmc.combined2 <- FindNeighbors(pbmc.combined2, dims = 1:9, verbose = FALSE)
pbmc.combined3 <- FindClusters(pbmc.combined2, resolution = 1, verbose = FALSE)
pbmc.combined3 <- RunUMAP(pbmc.combined3, dims = 1:9, verbose = FALSE) # 디멘션 건들면 모양 바뀜. 

plot2 <- DimPlot(pbmc.combined3, reduction = "umap", label=TRUE)
pbmc.combined3[["Sample_ID"]] <- pbmc.combined3$orig.ident

pbmc.combined3$celltype2 <- "HC"
pbmc.combined3$celltype2[(str_detect(names(pbmc.combined3$orig.ident), "GD"))] <- "S"
pbmc.combined3$celltype2[(str_detect(names(pbmc.combined3$orig.ident), "BD"))] <- "NS"
pbmc.combined3@meta.data$celltype2 <- factor(pbmc.combined3@meta.data$celltype2, levels= c("HC", "S", "NS"))

Idents(pbmc.combined3) <- pbmc.combined3$celltype2
plot3 <- DimPlot(pbmc.combined3, reduction = "umap", split.by = "celltype2", label=F)

saveRDS(pbmc.combined3, file = "D:/HP_rdata/GSE167363/6h_Sepsis_merging_seurat_data(cluster0_27).rds")
pbmc.combined3 <- readRDS(file = "D:/HP_rdata/GSE167363/6h_Sepsis_merging_seurat_data(cluster0_27).rds")

Idents(pbmc.combined3) <- pbmc.combined3$RNA_snn_res.1

### Cell type genes ###

ref1 <- HumanPrimaryCellAtlasData()
# ref2 <- NovershternHematopoieticData()
#ref3 <- DatabaseImmuneCellExpressionData()
#ref4 <- MonacoImmuneData()
pbmc.combined4 <- JoinLayers(pbmc.combined3)
DefaultAssay(pbmc.combined4) <- "RNA"
combine.all.sce <- as.SingleCellExperiment(pbmc.combined4)
pred.Result <- SingleR(test = combine.all.sce, ref = ref1, assay.type.test=1, labels = ref1$label.fine)
tab <- table(cluster = pbmc.combined4$seurat_clusters, label = pred.Result$labels)
pheatmap::pheatmap(log2(tab+1), cluster_cols = F)

# Labeling ---------------------------------------------------------------------
# 참조를 위해 이전 ID 클래스(클러스터 레이블)를 저장합니다.
Idents(pbmc.combined4)
Idents(pbmc.combined4) <- pbmc.combined4$seurat_clusters

pbmc.combined4 <- AddModuleScore(
  object = pbmc.combined4,
  features = feature4,
  ctrl = 5,
  name = 'Score(1)'
)

pbmc.combined4 <- AddModuleScore(
  object = pbmc.combined4,
  features = feature5,
  ctrl = 5,
  name = 'Score(2)'
)

boxplot(pbmc.combined4$`Score(1)1`)
mean(pbmc.combined4$`Score(2)1`)
mean(pbmc.combined4$`Score(1)1`)

max(pbmc.combined4$`Score(1)1`)
max(pbmc.combined4$`Score(1)2`)
max(pbmc.combined4$`Score(1)3`)
max(pbmc.combined4$`Score(2)1`)
pbmc.combined4$`Score(1)1`

FeaturePlot(pbmc.combined4, features = "Score(1)1")
FeaturePlot(pbmc.combined4, features = "Score(1)2")
FeaturePlot(pbmc.combined4, features = "Score(1)3")
FeaturePlot(pbmc.combined4, features = "Score(2)1")

pbmc.combined4$orig.ident

testpbmc <- subset(pbmc.combined4, ident = c("25", "13", "11", "9", "12", "21", "23", "15"))

ggplot(testpbmc@meta.data, aes(x=celltype2, fill=seurat_clusters)) + geom_bar(colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()

ggplot(testpbmc@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar(colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()


DimPlot(pbmc.combined4, split.by = "celltype2", label = T)

RidgePlot(pbmc.combined4, features = "Score(1)1", group.by = "seurat_clusters")

pbmc.combined4 <- RenameIdents(
  object = pbmc.combined4,
  `0` = "T cells",
  `1` = "NK cells",
  `2` = "T cells",
  `3` = "Monocytes",
  `4` = "T cells",
  `5` = "B cells",
  `6` = "T cells",
  `7` = "B cells",
  `8` = "B cells",
  `9` = "Macrophages",
  `10` = "Platelets",
  `11` = "Macrophages",
  `12` = "Macrophages",
  `13` = "Monocytes",
  `14` = "Platelets",
  `15` = "Neutrophils",
  `16` = "Erythoblasts",
  `17` = "Erythoblasts",
  `18` = "Monocytes",
  `19` = "Monocytes",
  `20` = "Undefined",
  `21` = "Macrophages",
  `22` = "Undefined",
  `23` = "Macrophages",
  `24` = "Monocytes",
  `25` = "Monocytes",
  `26` = "B cells",
  `27` = "Undefined"
)

pbmc.combined4$cell_annot <- Idents(pbmc.combined4)
Idents(pbmc.combined4) <- pbmc.combined4$cell_annot
plot4 <- DimPlot(pbmc.combined4, reduction = "umap", split.by = "celltype2", label=TRUE)
testpbmc <- subset(pbmc.combined4, ident = c("Macrophages"))
Idents(testpbmc) <- testpbmc$celltype2
testpbmc2 <- subset(testpbmc, ident = c("S", "NS"))
plot5 <- DimPlot(testpbmc2, reduction = "umap", split.by = "celltype2", label=TRUE)
Idents(testpbmc2) <- testpbmc2$seurat_clusters
FeaturePlot(testpbmc2, features = "CLEC4E", split.by = "seurat_clusters")

testpbmc2.markers <- FindAllMarkers(testpbmc2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(testpbmc2.markers, file = "D:/HP_rdata/GSE167363/6h_Macrophages_S_NS_markers.csv")
MoMa3.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

RidgePlot(MoMa5, features = target) + NoLegend()
DoHeatmap(MoMa3, features = target) + NoLegend()


ggplot(testpbmc@meta.data, aes(x=celltype2, fill=seurat_clusters)) + geom_bar(colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()


pbmc.combined4$cell_annot <- Idents(pbmc.combined4)
pbmc.combined4$orig.ident
cols <- colorRampPalette(c("#EF9A9A", "#03A9F4", "#AED581"))(11)
ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=cell_annot)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()
pbmc.combined4$orig.ident <- factor(pbmc.combined4$orig.ident, levels = c("HC1", "HC2", "GD1",))
VlnPlot(pbmc.combined4, feature = "HP", group.by = "orig.ident")

Idents(pbmc.combined4) <- pbmc.combined4$cell_annot
DimPlot(pbmc.combined4, split.by = "celltype2")

Idents(pbmc.combined4) <- pbmc.combined4$celltype2
DimPlot(pbmc.combined4, split.by = "celltype2")
FeaturePlot(pbmc.combined4, feature = c("HP", "FUT4", "CLEC4E", "IL1B"),  split.by = "celltype2")

MoMa1 <- subset(pbmc.combined4, ident = c("4", "13", "20", "22", "8"))
Idents(MoMa1) <- MoMa1$celltype2
MoMa2 <- subset(MoMa1, ident = c("S", "NS"))
Idents(MoMa2) <- MoMa2$celltype2
MoMa2$seurat_clusters

cols <- rev(colorRampPalette(c("#EF9A9A", "#03A9F4", "#AED581"))(9))
ggplot(MoMa2@meta.data, aes(x=celltype2, fill=seurat_clusters)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()

feature1 <- list(c("HP", "CLEC4E"))
feature2 <- list(c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9", "FUT10", "FUT11", "FUT3", "FUT3", "POFUT1", "POFUT2"))

cclcxcl <- read.csv(file = "d:/GO_term/cytokine_chemokine.csv")
cclcxcl2 <- c(cclcxcl$Cytokine, cclcxcl$Chemokine)
cclcxcl3 <- cclcxcl2[!duplicated(cclcxcl2)]
cclcxcl4 <- cclcxcl3[c(1,2,7,9,11,17, 22:67)]

feature3 <- list(cclcxcl4)
feature4 <- list(c("HP", "CLEC4E"), c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9", "FUT10", "FUT11", "FUT3", "FUT3", "POFUT1", "POFUT2"), 
                 cclcxcl4)

feature5 <- list(c("HP", "CLEC4E", "IL1B", "FUT4"))

write.csv(cclcxcl4, file = "d:/GO_term/cclcxcl.csv")

pbmc.combined5 <- AddModuleScore(
  object = pbmc.combined5,
  features = feature5,
  ctrl = 5,
  name = 'HP_CLEC4E'
)
head(x = pbmc.combined5[])

ggplot(pbmc.combined5@meta.data, aes(x=orig.ident, fill=Prognosis1)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()


pbmc.combined4$orig.ident <- factor(pbmc.combined4$orig.ident, levels = c("HC1", "HC2", "GD1", "GD2", "GD3", "BD1", "BD2"))

FeaturePlot(pbmc.combined5, feature = c("HP_CLEC4E1"),  split.by = "celltype2")
RidgePlot(pbmc.combined5, feature = c("Prognosis1", "Prognosis2", "Prognosis3"), group.by = "celltype2")

MoMa3 <- subset(pbmc.combined5, ident = c("25", "13", "11", "9", "12", "21", "23"))
FeaturePlot(MoMa3, feature = c("HP_CLEC4E1"),  split.by = "celltype2")
RidgePlot(MoMa3, feature = c("HP_CLEC4E1"), group.by = "celltype2")


Idents(pbmc.combined5) <- pbmc.combined5$seurat_clusters
MoMa3.markers <- FindAllMarkers(pbmc.combined5, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(MoMa3.markers, file = "D:/HP_rdata/GSE167363/S_NS_markers.csv")
MoMa3.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

RidgePlot(MoMa5, features = target) + NoLegend()
DoHeatmap(MoMa3, features = target) + NoLegend()

FeaturePlot(MoMa2, features = "IL1B")
DoHeatmap(MoMa2)
Idents(MoMa2)


ggplot(MoMa3@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar(colour = "black")

