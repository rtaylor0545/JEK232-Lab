## sepsis 환자 scRNAseq ##

library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

HC1 <- Read10X(data.dir = "D:/HP_rdata/HC1/")
HC1_so <- CreateSeuratObject(counts = HC1, project = "HC1")
HC1_so[["percent.mt"]] = PercentageFeatureSet(HC1_so, pattern = "^MT.")
HC2 <- Read10X(data.dir = "D:/HP_rdata/HC2/")
HC2_so <- CreateSeuratObject(counts = HC2, project = "HC2")
HC2_so[["percent.mt"]] = PercentageFeatureSet(HC2_so, pattern = "^MT.")

GD1 <- Read10X(data.dir = "D:/HP_rdata/GD1/")
GD1_so <- CreateSeuratObject(counts = GD1, project = "GD1")
GD1_so[["percent.mt"]] = PercentageFeatureSet(GD1_so, pattern = "^MT.")
GD2 <- Read10X(data.dir = "D:/HP_rdata/GD2/")
GD2_so <- CreateSeuratObject(counts = GD2, project = "GD2")
GD2_so[["percent.mt"]] = PercentageFeatureSet(GD2_so, pattern = "^MT.")

BD1 <- Read10X(data.dir = "D:/HP_rdata/BD1/")
BD1_so <- CreateSeuratObject(counts = BD1, project = "BD1")
BD1_so[["percent.mt"]] = PercentageFeatureSet(BD1_so, pattern = "^MT.")
BD2 <- Read10X(data.dir = "D:/HP_rdata/BD2/")
BD2_so <- CreateSeuratObject(counts = BD2, project = "BD2")
BD2_so[["percent.mt"]] = PercentageFeatureSet(BD2_so, pattern = "^MT.")

test <- readRDS(file = "d:/HP_rdata/RF230246_Seurat.rds")

# Merging all samples
pbmc.combined <- merge(HC1_so, y = c(HC2_so, GD1_so, GD2_so, BD1_so, BD2_so), 
                       add.cell.ids = c("HC1", "HC2", "GD1", "GD2", "BD1", "BD2"), 
                       project = "Sepsis")

table(pbmc.combined$orig.ident)

#Scatter 플랏
plot0 <- VlnPlot(pbmc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width = 12, repr.plot.height = 6)
plot1 + plot2

pbmc.combined <- subset(pbmc.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
pbmc.combined <- NormalizeData(pbmc.combined, verbose = FALSE)
pbmc.combined <- FindVariableFeatures(pbmc.combined, verbose = FALSE)
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)

saveRDS(pbmc.combined, file = "D:/HP_rdata/Sepsis_merging_seurat_data(prototype).rds")
pbmc.combined <- readRDS(file = "D:/HP_rdata/Sepsis_merging_seurat_data(prototype).rds")
write.csv(pbmc.combined@active.ident, file = "d:/cell_population.csv")

plot <- ElbowPlot(pbmc.combined)
options(repr.plot.width = 6, repr.plot.height = 6)
plot

pbmc.combined2 <- FindNeighbors(pbmc.combined, dims = 1:9, verbose = FALSE)
pbmc.combined3 <- FindClusters(pbmc.combined2, resolution = 0.5, verbose = FALSE)
pbmc.combined3 <- RunUMAP(pbmc.combined3, dims = 1:9, verbose = FALSE) # 디멘션 건들면 모양 바뀜. 
pbmc.combined3 <- RunTSNE(pbmc.combined3, dims = 1:9, verbose = FALSE)

plot2 <- DimPlot(pbmc.combined3, reduction = "umap", label=TRUE)
plot3 <- DimPlot(pbmc.combined3, reduction = "tsne", label=TRUE)

pbmc.combined3[["seurat_cluster"]] <- pbmc.combined3@active.ident
pbmc.combined3[["Sample_ID"]] <- pbmc.combined$orig.ident

saveRDS(pbmc.combined3, file = "D:/HP_rdata/Sepsis_merging_seurat_data(cluster0_15).rds")
pbmc.combined3 <- readRDS(file = "D:/HP_rdata/Sepsis_merging_seurat_data(cluster0_15).rds")

table(pbmc.combined3$orig.ident)

# Cell type dot plot 

# 0 - NK
# 1 - CD4T
# 2 - CD8T
# 3 - B cell
# 4 - CD4T
# 5 - mono
# 6 - NK
# 7 - Macrophage
# 8 - Macrophage
# 9 - CD16+ Mono 
# 10 - CD16+ Mono 
# 11 - CD14 mono?? Macrophage ?? ????
# 12 - nk
# 13 - Macrophage
# 14 - B cell
# 15 - Macrophage

CD4T <-  c("CD3D", "CD3E", "CD4", "IL2RA", "LCK")
Bcells <- c("CD19", "CD79A", "MS4A1", "BANK1", "HLA-DRA")  # 3
NK <- c("NCAM1", "KIR2DL1", "KIR2DL2", "KIR2DL3", "NKG7", "PRF1", "GZMB")  # 0
CD8T <- c("CD3D", "CD3E", "CD8A", "CD8B", "GZMK")
CD14Mono <- c("CD14", "LYZ", "S100A8", "S100A9", "FPR1")
CD16Mono <- c("FCGR3A", "CD14", "S100A8", "S100A9", "LYZ")
NP <- c("FCGR3B", "ELANE", "MPO", "PRTN3", "CXCR1")
MP <- c("CD14", "FCGR3A", "MIR155HG", "CCL2", "IL1A", "IL1B", "IL6", "NOS2", "TLR2", "TLR4", "CD80", "CD86", "IL10", "IL4", "IL13", "CD115", "CD206", "PPARG", "ARG1", "CD163","CD68", "CD84", "CSF1R", "AIF1", "C1QA")
Effector <- c("CXCR5", "IL21", "IL17")
total <- c(CD4T, Bcells, NK, CD8T, CD14Mono, CD16Mono, NP, MP)
total2 <- total[!duplicated(total)]
p2 <- DotPlot(pbmc.combined4, features = Effector)

library(dplyr)
library(stringr)

pbmc.combined3$celltype <- "HC"
pbmc.combined3$celltype[!(str_detect(names(pbmc.combined$orig.ident), "HC"))] <- "SS"
pbmc.combined3$celltype

Idents(pbmc.combined3) <- "celltype"
plot5 <- DimPlot(pbmc.combined3, reduction = "umap", split.by = "celltype")

### Cell type genes ###
library(celldex)
library(SingleR)
library(dplyr)

ref1 <- HumanPrimaryCellAtlasData()
# ref2 <- NovershternHematopoieticData()
ref3 <- DatabaseImmuneCellExpressionData()
ref4 <- MonacoImmuneData()
pbmc.combined4 <- JoinLayers(pbmc.combined3)
DefaultAssay(pbmc.combined4) <- "RNA"
combine.all.sce <- as.SingleCellExperiment(pbmc.combined4)
pred.Result <- SingleR(test = combine.all.sce, ref = ref4, assay.type.test=1, labels = ref4$label.fine)
tab <- table(cluster = pbmc.combined4$seurat_clusters, label = pred.Result$labels)
pheatmap::pheatmap(log10(tab+10), cluster_cols = F)

pbmc.combined4$Sample_ID

# Labeling ---------------------------------------------------------------------
# 참조를 위해 이전 ID 클래스(클러스터 레이블)를 저장합니다.
Idents(pbmc.combined4) <- "seurat_cluster"

# 레이블 변경하기
pbmc.combined4 <- RenameIdents(
  object = pbmc.combined4,
  `0` = "NK cells",
  `1` = "CD4+ T cells",
  `2` = "CD4+ T cells",
  `3` = "B cells",
  `4` = "CD4+ T cells",
  `5` = "Monocytes",
  `6` = "CD8+ T cells",
  `7` = "Macrophages-1",
  `8` = "Macrophages-2", 
  `9` = "Monocytes", 
  `10` = "Monocytes", 
  `11` = "Macrophages-1",
  `12` = "NK cells",
  `13` = "Macrophages-2",
  `14` = "B cells",
  `15` = "Macrophages-2"
)

p <- DimPlot(pbmc.combined4, label = T, reduction = "umap", split.by = "celltype")
options(repr.plot.width = 7, repr.plot.height = 7)

pbmc.combined4@active.ident <- 
  factor(pbmc.combined4@active.ident, 
         levels = rev(c("NK cells", "CD4+ T cells", "CD8+ T cells", "Macrophages-1", "Macrophages-2", "Monocytes", "B cells")))

target <- c("GZMB", "IL7R", "CD247", "HAVCR2", "SELL", "LYZ", "HLA-DPA1", "CD79A")

p3 <- DotPlot(pbmc.combined4, features = target)

Idents(pbmc.combined4) <- "celltype"
HC <- subset(pbmc.combined4, ident = c("HC"))
SS <- subset(pbmc.combined4, ident = c("SS"))

Idents(HC) <- "seurat_cluster"
Idents(SS) <- "seurat_cluster"

SS <- RenameIdents(
  object = SS,
  `0` = "NK cells",
  `1` = "CD4+ T cells",
  `2` = "CD4+ T cells",
  `3` = "B cells",
  `4` = "CD4+ T cells",
  `5` = "Monocytes",
  `6` = "CD8+ T cells",
  `7` = "Macrophages-1",
  `8` = "Macrophages-2", 
  `9` = "Monocytes", 
  `10` = "Monocytes", 
  `11` = "Macrophages-1",
  `12` = "NK cells",
  `13` = "Macrophages-2",
  `14` = "B cells",
  `15` = "Macrophages-2"
)

Idents(HC) <- "celltype"

IF <- c("IL12B", "CXCL11", "CXCL10", "CXCL2", "CXCL3", "CXCL5", "CXCL1", "CXCL6", 
        "CXCL8", "CCL20", "CCL18", "CCL2", "CCL3", "CCL4", "IL1B", "TNF")
FTF2 <- c(paste0(rep("FUT", 9), seq(3, 11, 1)), "HP")
VlnPlot(HC, features = FTF2, split.by = "celltype", pt.size = 0)
VlnPlot(SS, features = FTF2, split.by = "celltype", pt.size = 0)
VlnPlot(HC, features = IF, pt.size = 0)
VlnPlot(SS, features = IF, pt.size = 0)

saveRDS(pbmc.combined4, file = "Celltype_Sepsis_seurat_data(p4).rds")
pbmc.combined4 <- readRDS(file = "d:/HP_rdata/Celltype_Sepsis_seurat_data(p4).rds")

p3 <- DotPlot(pbmc.combined4, features = "TFEB")
p4 <- VlnPlot(pbmc.combined4, features = "TFEB", split.by = "celltype", pt.size = 0)
pbmc.combined4$seurat_cluster

pbmc.combined4$Cell_anot <- Idents(pbmc.combined4)
Idents(pbmc.combined4) <- pbmc.combined4$orig.ident
pbmc.combined5 <- subset(pbmc.combined4, ident = c("HC1", "HC2"))
pbmc.combined6 <- subset(pbmc.combined4, ident = c("BD1", "BD2"))
pbmc.combined7 <- subset(pbmc.combined4, ident = c("BD1"))
pbmc.combined8 <- subset(pbmc.combined4, ident = c("BD2"))
table(pbmc.combined5$Cell_anot)
table(pbmc.combined6$Cell_anot)
table(pbmc.combined7$Cell_anot)
table(pbmc.combined8$Cell_anot)
Idents(pbmc.combined5)
Idents(pbmc.combined5) <- pbmc.combined5$Cell_anot
VlnPlot(pbmc.combined5, features = c("HP", "FUT4", "CLEC4E", "IL1B"), pt.size = 0, ncol = 4)
HC <- subset(pbmc.combined5, ident = c("Monocytes"))
Ma1 <- subset(pbmc.combined5, ident = c("Macrophages-1", "Macrophages-2"))
table(Mo1$Cell_anot)

MoMa <- readRDS(file = "d:/HP_rdata/S_NS_Mph(240321Final).rds")
MoMa <- subset(pbmc.combined5, ident = c("Macrophages-1", "Macrophages-2"))
MoMa.markers <- FindAllMarkers(MoMa, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5)
write.csv(MoMa.markers, file = "d:/HP_rdata/MoMA_FC_markers(0.05_0.5).csv")


pbmc.combined4$seurat_clusters

table(pbmc.combined4$orig.ident)

pbmc.combined4$Cell_annot <- Idents(pbmc.combined4)

md <- pbmc.combined4@meta.data %>% as.data.table
head(AverageExpression(object = pbmc.combined4, group.by = c('ident', 'groups'))$RNA)

pbmc.markers <- FindAllMarkers(pbmc.combined4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
write.csv(pbmc.markers, file = "d:/HP_rdata/celltype_markers.csv")

#basic plot of clusters by replicate
ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=Cell_annot)) + geom_bar()

#set custom order to dataset w 30 clusters
my_levels <- c(0, 12, 1, 2, 4, 6, 7, 11, 8, 13, 15, 5, 9, 10, 3, 14)
factor(pbmc.combined4@meta.data$RNA_snn_res.0.5, levels= my_levels)
pbmc.combined4@meta.data$RNA_snn_res.0.5 <- factor(pbmc.combined4@meta.data$RNA_snn_res.0.5, levels= my_levels)
ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=RNA_snn_res.0.5)) + geom_bar(position = "fill")
ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=RNA_snn_res.0.5)) + geom_bar(position = "fill") +
scale_fill_manual(values=cols)

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Blues")[7]))(16))

#plot as proportion or percentage of cluster
ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=RNA_snn_res.0.5)) + geom_bar(position = "fill") +
  scale_fill_manual(values=cols)

display.brewer.all()

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greens")[7], 
                                brewer.pal(n = 9, name = "Blues")[7]))(7))

cols <- rev(colorRampPalette(c("#CA7797", "#97CA77", "#7797CA"))(7))
cols <- c("#EF9A9A", "#F48FB1", "#CE93D8", "#B39DDB", "#9FA8DA", "#90CAF9", "#80DEEA")
cols <- c("#C62828", "#AD1457", "#6A1B9A", "#4527A0", "#283593", "#1565C0", "#00838F")
cols <- rev(colorRampPalette(c("#FFEBEE", "#E3F2FD", "#F1F8E9"))(7))
cols <- rev(colorRampPalette(c("#EF9A9A", "#03A9F4", "#AED581"))(7))

my_levels <- c("NK cells", "CD4+ T cells", "CD8+ T cells", "B cells", "Monocytes", "Macrophages-1", "Macrophages-2")
factor(pbmc.combined4@meta.data$RNA_snn_res.0.5, levels= my_levels)
pbmc.combined4@meta.data$Cell_annot <- factor(pbmc.combined4@meta.data$RNA_snn_res.0.5, levels= my_levels)

ggplot(pbmc.combined4@meta.data, aes(x=orig.ident, fill=Cell_annot)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols) +coord_flip()

Idents(pbmc.combined4)
MoMa1 <- subset(pbmc.combined4, ident = c("Monocytes", "Macrophages-1", "Macrophages-2"))
MoMa1$celltype
MoMa1$Group <- "Monocytes"
MoMa1$Group[!(str_detect(MoMa1$Cell_annot, "Mono"))] <- "Macrophages"
factor(MoMa1$Group)

MoMa1@meta.data$Group <- factor(MoMa1@meta.data$Group, levels= c("Monocytes", "Macrophages"))
MoMa1@meta.data$celltype <- factor(MoMa1@meta.data$celltype, levels= c("HC", "SS"))
saveRDS(MoMa1, file = "d:/HP_rdata/Monocytes_Macrophages_seurat.rds")
MoMa1 <- readRDS(file = "d:/HP_rdata/Monocytes_Macrophages_seurat.rds")
rm(p6)
ggplot(MoMa1@meta.data, aes(x=celltype, fill=Cell_annot, label = scales::percent(Cell_annot))) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols)

HC_All <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$celltype == "HC")])
HC_Mo <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$Cell_annot == "Monocytes") & (MoMa1@meta.data$celltype == "HC")])
HC_Ma1 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages-1" &  (MoMa1@meta.data$celltype == "HC")])
HC_Ma2 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages-2" &  (MoMa1@meta.data$celltype == "HC")])

HC_Mo / HC_All * 100
HC_Ma1 / HC_All * 100
HC_Ma2 / HC_All * 100

SS_All <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$celltype == "SS")])
SS_Mo <- length(MoMa1@meta.data$Cell_annot[(MoMa1@meta.data$Cell_annot == "Monocytes") & (MoMa1@meta.data$celltype == "SS")])
SS_Ma1 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages-1" &  (MoMa1@meta.data$celltype == "SS")])
SS_Ma2 <- length(MoMa1@meta.data$Cell_annot[MoMa1@meta.data$Cell_annot == "Macrophages-2" &  (MoMa1@meta.data$celltype == "SS")])

SS_Mo / SS_All * 100
SS_Ma1 / SS_All * 100
SS_Ma2 / SS_All * 100

Idents(MoMa1)
MoMa1.markers <- FindAllMarkers(MoMa1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
write.csv(MoMa1.markers, file = "d:/HP_rdata/MoMa_celltype_markers.csv")
MoMa2 <- subset(MoMa1, ident = c("Macrophages-1", "Macrophages-2"))
MoMa2.markers <- FindAllMarkers(MoMa2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
MoMa2.markers2 <- FindAllMarkers(MoMa2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
write.csv(MoMa2.markers, file = "d:/HP_rdata/Mph1_Mph2_markers.csv")
write.csv(MoMa2.markers2, file = "d:/HP_rdata/Mph1_Mph2_markers2.csv")
VlnPlot(MoMa2, features = c("HP", "FUT4", "CLEC4E", "IL1B"), pt.size = 0, ncol = 4)
VlnPlot(MoMa2, features = c("CXCL2", "CXCL3", "CXCL8", "CCL3"), pt.size = 0, ncol = 4)

p <- DimPlot(MoMa2, label = T, reduction = "umap")
options(repr.plot.width = 7, repr.plot.height = 7)
p
MoMa2.markers <- FindAllMarkers(MoMa2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
write.csv(MoMa2.markers, file = "d:/HP_rdata/Ma2_celltype_markers.csv")

MoMa2.markers <- read.csv(file = "d:/HP_rdata/Ma2_celltype_markers.csv", row.names = 1)

MoMa2.markers %>%
  group_by(cluster) %>%
  slice_max(n = 7, order_by = avg_log2FC) -> top10

test <- c("HP", "FUT4", "CLEC4E", "IL1B")

DoHeatmap(MoMa2, features = test) + NoLegend()

FeaturePlot(MoMa2, features = "IL1B")
DoHeatmap(MoMa2)
Idents(MoMa2)
