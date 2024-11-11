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


HC1 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/HC1/")
HC1_so <- CreateSeuratObject(counts = HC1, project = "HC1")
HC1_so[["percent.mt"]] = PercentageFeatureSet(HC1_so, pattern = "^MT.")
HC2 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/HC2/")
HC2_so <- CreateSeuratObject(counts = HC2, project = "HC2")
HC2_so[["percent.mt"]] = PercentageFeatureSet(HC2_so, pattern = "^MT.")

GD1 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/SS1/")
GD1_so <- CreateSeuratObject(counts = GD1, project = "GD1")
GD1_so[["percent.mt"]] = PercentageFeatureSet(GD1_so, pattern = "^MT.")
GD2 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/SS2/")
GD2_so <- CreateSeuratObject(counts = GD2, project = "GD2")
GD2_so[["percent.mt"]] = PercentageFeatureSet(GD2_so, pattern = "^MT.")
GD3 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/SS3/")
GD3_so <- CreateSeuratObject(counts = GD3, project = "GD3")
GD3_so[["percent.mt"]] = PercentageFeatureSet(GD3_so, pattern = "^MT.")

BD1 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/NS1/")
BD1_so <- CreateSeuratObject(counts = BD1, project = "BD1")
BD1_so[["percent.mt"]] = PercentageFeatureSet(BD1_so, pattern = "^MT.")
BD2 <- Read10X(data.dir = "D:/HP_rdata/GSE167363/GSE167363_RAW/NS2/")
BD2_so <- CreateSeuratObject(counts = BD2, project = "BD2")
BD2_so[["percent.mt"]] = PercentageFeatureSet(BD2_so, pattern = "^MT.")

# Merging all samples
pbmc.combinedd1 <- merge(HC1_so, y = c(HC2_so, GD1_so, GD3_so, BD1_so, BD2_so), 
                       add.cell.ids = c("HC1", "HC2", "GD1", "GD3", "BD1", "BD2"), 
                       project = "GSE167363")

table(pbmc.combinedd1$orig.ident)

#Scatter 플랏
plot0 <- VlnPlot(pbmc.combinedd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc.combinedd1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.combinedd1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width = 12, repr.plot.height = 6)
plot0 
plot1 + plot2

rm(plot2)

pbmc.combined2 <- subset(pbmc.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mt < 25)
pbmc.combined2 <- NormalizeData(pbmc.combined2, verbose = FALSE)
pbmc.combined2 <- FindVariableFeatures(pbmc.combined2, verbose = FALSE)
pbmc.combined2 <- ScaleData(pbmc.combined2, verbose = FALSE)
pbmc.combined2 <- RunPCA(pbmc.combined2, verbose = FALSE)

saveRDS(pbmc.combined, file = "D:/HP_rdata/GSE167363/Sepsis_merging_seurat_data(prototype).rds")
saveRDS(pbmc.combined2, file = "D:/HP_rdata/GSE167363/Sepsis_merging_seurat_data(filter).rds")
write.csv(pbmc.combined@active.ident, file = "d:/HP_rdata/GSE167363/cell_population(proto).csv")
write.csv(pbmc.combined2@active.ident, file = "d:/HP_rdata/GSE167363/cell_population(filter).csv")

############# Analysis ###############
plot <- ElbowPlot(pbmc.combined2)
options(repr.plot.width = 6, repr.plot.height = 6)
plot

rm(plot)

pbmc.combined2 <- FindNeighbors(pbmc.combined2, dims = 1:9, verbose = FALSE)
pbmc.combined3 <- FindClusters(pbmc.combined2, resolution = 1, verbose = FALSE)
pbmc.combined3 <- RunUMAP(pbmc.combined3, dims = 1:9, verbose = FALSE) # 디멘션 건들면 모양 바뀜. 
#pbmc.combined3 <- RunTSNE(pbmc.combined3, dims = 1:9, verbose = FALSE)

plot2 <- DimPlot(pbmc.combined3, reduction = "umap", label=TRUE)
#plot3 <- DimPlot(pbmc.combined3, reduction = "tsne", label=TRUE)

pbmc.combined3[["Sample_ID"]] <- pbmc.combined3$orig.ident
pbmc.combined3$celltype <- "HC"
pbmc.combined3$celltype[!(str_detect(names(pbmc.combined3$orig.ident), "HC"))] <- "SS"

saveRDS(pbmc.combined3, file = "D:/HP_rdata/GSE167363/Sepsis_merging_seurat_data(cluster0_30).rds")
pbmc.combined3 <- readRDS(file = "D:/HP_rdata/GSE167363/Sepsis_merging_seurat_data(cluster0_30).rds")

Idents(pbmc.combined3) <- pbmc.combined3$celltype
plot3 <- DimPlot(pbmc.combined3, label = F, reduction = "umap", split.by = 'celltype')

Idents(pbmc.combined3) <- pbmc.combined3$RNA_snn_res.0.5

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
pheatmap::pheatmap(log10(tab+10), cluster_cols = F)
pheatmap::pheatmap(log2(tab+1), cluster_cols = F)

# Cell type dot plot 
Idents(pbmc.combined4) <- "celltype"


# Labeling ---------------------------------------------------------------------
# 참조를 위해 이전 ID 클래스(클러스터 레이블)를 저장합니다.
pbmc.combined4 <- RenameIdents(
  object = pbmc.combined4,
  `24` = "Macrophages",
  `8` = "Macrophages",
  `6` = "Macrophages",
  `17` = "Macrophages",
  `19` = "Macrophages",
  `18` = "Monocytes",
  `25` = "Monocytes",
  `21` = "Monocytes",
  `4` = "Monocytes",
  `23` = "Monocytes",
  `0` = "T cells",
  `2` = "T cells",
  `13` = "T cells",
  `12` = "T cells",
  `9` = "NK cells",
  `5` = "NK cells",
  `10` = "NK cells",
  `1` = "B cells",
  `14` = "B cells",
  `3` = "B cells",
  `7` = "B cells",
  `16` = "Platelets",
  `20` = "Platelets",
  `11` = "Platelets",
  `22` = "Undefined",
  `26` = "Erythroblast",
  `27` = "Erythroblast",
  `15` = "Neutrophil",
  `28` = "Undefined",
  `29` = "Platelets",
  `30` = "Undefined"
)
pbmc.combined4$cellanot <- Idents(pbmc.combined4)

p <- DimPlot(pbmc.combined4, label = T, reduction = "umap", split.by = "celltype")

pbmc.combined4$celltype2 <- "HC"
pbmc.combined4$celltype2[(str_detect(names(pbmc.combined4$orig.ident), "GD"))] <- "S"
pbmc.combined4$celltype2[(str_detect(names(pbmc.combined4$orig.ident), "BD"))] <- "NS"
pbmc.combined4@meta.data$celltype2 <- factor(pbmc.combined4@meta.data$celltype2, levels= c("HC", "S", "NS"))

p2 <- DimPlot(pbmc.combined4, label = T, reduction = "umap", split.by = "celltype2")
saveRDS(pbmc.combined4, file = "D:/HP_rdata/GSE167363/Sepsis_merging_seurat_data(cluster0_30_Annotation_Notused).rds")

### 여기서부터 시작하시면 됨. ###--------------------------------------------------------------------------------------------------

pbmc.combined4 <- readRDS(file = "D:/HP_rdata/GSE167363/Sepsis_merging_seurat_data(cluster0_30_Annotation_Notused).rds")
Idents(pbmc.combined4) <- pbmc.combined4$cellanot
pbmc.combined4@meta.data$celltype2 <- factor(pbmc.combined4@meta.data$celltype2, levels= c("HC", "S", "NS"))
plot2 <- DimPlot(pbmc.combined4, reduction = "umap", label=T, split.by = "celltype2")

pbmc.combined4$celltype
pbmc.combined4$cellanot
pbmc.combined4$seurat_clusters
pbmc.combined4$cellanot

pbmc.combined4$orig.ident[(str_detect(names(pbmc.combined4$orig.ident), "GD3"))] <- "GD3"
feature5 <- list(c("HP", "CLEC4E", "IL1B", "FUT4"))

pbmc.combined4 <- AddModuleScore(
  object = pbmc.combined4,
  features = feature5,
  ctrl = 5,
  name = 'Score(5)'
)
pbmc.combined4$`Score(5)1`
pbmc.combined4$celltype3 <- "X"
pbmc.combined4$celltype3[pbmc.combined4$`Score(5)1` > 0] <- "O"
pbmc.combined4$celltype3 <- factor(pbmc.combined4@meta.data$celltype3, levels= c("X", "O"))
cols1 <- rev(colorRampPalette(c("#F40303", "#FFFFFF"))(7))
ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=celltype3)) + geom_bar(position = "fill") +
  scale_fill_manual(values=cols1)

Idents(pbmc.combined4) <- pbmc.combined4$celltype3
pbmc.combined5 <- subset(pbmc.combined4, ident = "O")

ggplot(pbmc.combined5@meta.data, aes(x=celltype2, fill=celltype3)) + geom_bar(color = "black") +
  scale_fill_manual(values=cols1)

table(pbmc.combined5$orig.ident)

DimPlot(pbmc.combined4, reduction = "umap", label=T, split.by = "celltype2")
target <- c("IL1B", "FUT4", "HP", "CLEC4E")
DotPlot(pbmc.combined4, features = target)

cols1 <- rev(colorRampPalette(c("#F40303", "#FFFFFF", "#0006FF"))(31))

ggplot(pbmc.combined4@meta.data, aes(x=celltype2, fill=RNA_snn_res.1)) + geom_bar(position = "fill", linewidth = 5, linetype = 3)
Idents(pbmc.combined4) <- pbmc.combined4$RNA_snn_res.1
Ma24 <- subset(pbmc.combined4, ident = "24")
ggplot(Ma24@meta.data, aes(x=celltype2, fill=RNA_snn_res.1)) + geom_bar(linewidth = 5, linetype = 3)

pbmc.combined4$celltype3 <- "X"
pbmc.combined4$celltype3[pbmc.combined4$`Score(5)1` > 0] <- "O"
Idents(pbmc.combined4) <- pbmc.combined4$celltype3
pbmc.combined5 <- subset(pbmc.combined4, ident = "O")
Idents(pbmc.combined5) <- pbmc.combined5$celltype2
VlnPlot(pbmc.combined5, features = c("score2"))


Ma1 <- subset(pbmc.combined4, ident = "Macrophages")
pbmc.combined4$score2 <- pbmc.combined4$`Score(5)1`
pbmc.combined4$score2[pbmc.combined4$score2 < 0] <- 0
VlnPlot(pbmc.combined4, features = c("score2"))


library(tidyverse)
library(patchwork)
library(viridis)

p1 <- FeaturePlot(Ma1, "Score(5)1", order = TRUE, split.by = "celltype2") & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 2.5))
p2 <- FeaturePlot(pbmc.combined4, "Score(5)1", order = TRUE) & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 10))
wrap_plots(p1, p2)

FeaturePlot(pbmc.combined4, features = "Score(5)1", split.by = "celltype2") & scale_color_gradientn(colors = plasma(n = 10, direction = -1), limits = c(0, 10))
VlnPlot(pbmc.combined4, features = c("Score(5)1"))

table(pbmc.combined4$orig.ident[pbmc.combined4$`Score(5)1` > 0.5])

Idents(pbmc.combined4) <- pbmc.combined4$cellanot
Ma1 <- subset(pbmc.combined4, ident = "Macrophages")
table(Ma1$orig.ident)
table(Ma1$orig.ident[Ma1$`Score(2)1` > 0.5])
ta1 <- as.data.frame(NA)
ta1$BD1 <- 56/266
ta1$BD2 <- 20/331
ta1$GD1 <- 132/1452
ta1$GD2 <- 23/121
ta1$GD3 <- 55/1770

head(x = pbmc.combined4[])
FeaturePlot(Ma1, features = "Score(2)1", split.by = "orig.ident", pt.size = 0.5, cols = c(colorRampPalette(c("gray", "red"))(5)))

table(pbmc.combined4$orig.ident)
table(pbmc.combined4$cellanot)

cellnum1 <- table(pbmc.combined4$orig.ident[pbmc.combined4$cellanot == "Macrophages"])
cellnum1 <- c(cellnum1, HC1 = 0)
cellnum1 <- cellnum1[rev(order(names(cellnum1)))]

names(cellnum1) <- factor(names(cellnum1), levels = c("HC1", "HC2", "GD1", "GD2", "GD3", "BD1", "BD2"))
cellnum2 <- table(pbmc.combined4$orig.ident[pbmc.combined4$cellanot == "Monocytes"])
cellnum2 <- c(BD1 = 0, cellnum2)
names(cellnum2) <- factor(names(cellnum2), levels = c("HC1", "HC2", "GD1", "GD2", "GD3", "BD1", "BD2"))
cellnum3 <- table(pbmc.combined4$orig.ident[pbmc.combined4$`Score(2)1` > 0.2])
names(cellnum3) <- factor(names(cellnum3), levels = c("HC1", "HC2", "GD1", "GD2", "GD3", "BD1", "BD2"))
cellnum4 <- rbind(cellnum1, cellnum2, cellnum3)

testpbmc <- subset(pbmc.combined4, ident = c("Macrophages"))
DimPlot(testpbmc, reduction = "umap", label=TRUE, feature = "")
cols <- rev(colorRampPalette(c("#EF9A9A", "#03A9F4", "#AED581"))(9))

cols = c(colorRampPalette(c("red", "red"))(5))

ggplot(testpbmc@meta.data, aes(x=celltype2, fill=seurat_clusters)) + geom_bar(colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()

ggplot(testpbmc@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar(colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()





my_levels <- c("Macrophages", "Monocytes", "T cells", "B cells", "NK cells", "Platelets", "Undefined", "Erythroblast", "Neutrophil")
factor(pbmc.combined4@meta.data$Cell_annot, levels= my_levels)
pbmc.combined4@meta.data$Cell_annot <- factor(pbmc.combined4@meta.data$Cell_annot, levels= my_levels)

ggplot(pbmc.combined4@meta.data, aes(x=celltype2, fill=Cell_annot)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()

MoMa1 <- subset(pbmc.combined4, ident = c("Monocytes", "Macrophages"))
Idents(MoMa1) <- MoMa1$celltype2
factor(MoMa1$celltype2)

HC_All <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$celltype == "HC")])
HC_Mo <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$Cell_annot == "Monocytes") & (MoMa1@meta.data$celltype == "HC")])
HC_Ma1 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages" &  (MoMa1@meta.data$celltype == "HC")])

HC_Mo / HC_All * 100
HC_Ma1 / HC_All * 100

NS_All <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$celltype2 == "NS")])
NS_Mo <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$Cell_annot == "Monocytes") & (MoMa1@meta.data$celltype2 == "NS")])
NS_Ma1 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages" &  (MoMa1@meta.data$celltype2 == "NS")])

NS_Mo / NS_All * 100
NS_Ma1 / NS_All * 100

S_All <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$celltype2 == "S")])
S_Mo <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$Cell_annot == "Monocytes") & (MoMa1@meta.data$celltype2 == "S")])
S_Ma1 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages" &  (MoMa1@meta.data$celltype2 == "S")])

S_Mo / S_All * 100
S_Ma1 / S_All * 100

MoMa2 <- subset(pbmc.combined4, ident = c("Macrophages"))
Idents(MoMa2) <- MoMa2$RNA_snn_res.1
MoMa3 <- subset(MoMa2, ident = c("S", "NS"))
MoMa3.markers <- FindAllMarkers(MoMa3, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(MoMa3.markers, file = "D:/HP_rdata/GSE167363/S_NS_markers.csv")
MoMa3.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

RidgePlot(MoMa5, features = target) + NoLegend()
DoHeatmap(MoMa3, features = target) + NoLegend()

FeaturePlot(MoMa2, features = "IL1B")
DoHeatmap(MoMa2)
Idents(MoMa2)

MoMa4 <- subset(MoMa2, ident = c("17", "19", "8", "6"))
MoMa5 <- subset(MoMa4, subset = HP > 0)
MoMa4$IL1B <- MoMa4@assays$RNA$scale.data["IL1B",]
MoMa4$HP <- MoMa4@assays$RNA$scale.data["HP",]
MoMa4$FUT4 <- MoMa4@assays$RNA$scale.data[,]
MoMa4$CLEC4E <- MoMa4@assays$RNA$scale.data["CLEC4E",]

Idents(MoMa5) <- MoMa5$celltype2

dim(MoMa4@assays$RNA$scale.data)

VlnPlot(MoMa3, features = c("HP", "CLEC4E", "FUT4", "IL1B"))
VlnPlot(MoMa5, features = c("HP"))

corr <- FindVariableFeatures(MoMa4, selection.method = "vst", nfeatures = 10000)
sd <- corr@assays$RNA$scale.data
cl8 <- names(MoMa4$seurat_clusters)[MoMa4$seurat_clusters == 8]
cl6 <- names(MoMa4$seurat_clusters)[MoMa4$seurat_clusters == 6]
cl17 <- names(MoMa4$seurat_clusters)[MoMa4$seurat_clusters == 17]
cl19 <- names(MoMa4$seurat_clusters)[MoMa4$seurat_clusters == 19]

Idents(pbmc.combined4) <- pbmc.combined4$seurat_clusters
MoMa1 <- subset(pbmc.combined4, ident = c("17", "19", "8", "6"))
Idents(MoMa1) <- MoMa1$celltype2
MoMa2 <- subset(MoMa1, ident = c("S", "NS"))
Idents(MoMa2) <- MoMa2$seurat_clusters
MoMa2.markers <- FindAllMarkers(MoMa2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(MoMa2.markers, file = "D:/HP_rdata/GSE167363/17-19-8-6_cluster_markers.csv")
