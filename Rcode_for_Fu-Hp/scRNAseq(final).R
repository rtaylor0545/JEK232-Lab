### Set 26 Un vs SP-Hp treated -------------------------------------------------

library(ggplot2)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(patchwork)
library(TENxPBMCData)
library(scDblFinder)
library(future)
library(pheatmap)

# Untreated samples

counts = read.table("d:/scRNAseq/HC/HC25.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
Untreated25 = CreateSeuratObject(counts = t(counts), project = "Untreated")
Untreated25[["percent.mt"]] = PercentageFeatureSet(Untreated25, pattern = "^MT.")

counts = read.table("d:/scRNAseq/HC/HC26.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
Untreated26 = CreateSeuratObject(counts = t(counts), project = "Untreated")
Untreated26[["percent.mt"]] = PercentageFeatureSet(Untreated26, pattern = "^MT.")

counts = read.table("d:/scRNAseq/HC/HC27.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
Untreated27 = CreateSeuratObject(counts = t(counts), project = "Untreated")
Untreated27[["percent.mt"]] = PercentageFeatureSet(Untreated27, pattern = "^MT.")

library(RColorBrewer)

display.brewer.all()
# SP-Hp samples

counts = read.table("d:/scRNAseq/SS/SS25.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
SP_Hp25 = CreateSeuratObject(counts = t(counts), project = "SP-Hp treated")
SP_Hp25[["percent.mt"]] = PercentageFeatureSet(SP_Hp25, pattern = "^MT.")

counts = read.table("d:/scRNAseq/SS/SS26.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
SP_Hp26 = CreateSeuratObject(counts = t(counts), project = "SP-Hp treated")
SP_Hp26[["percent.mt"]] = PercentageFeatureSet(SP_Hp26, pattern = "^MT.")

counts = read.table("d:/scRNAseq/SS/SS27.csv", skip = 5, sep = ",", header = TRUE, row.names = 1)
SP_Hp27 = CreateSeuratObject(counts = t(counts), project = "SP-Hp treated")
SP_Hp27[["percent.mt"]] = PercentageFeatureSet(SP_Hp27, pattern = "^MT.")


# Merging all samples

pbmc.combined <- merge(Untreated25, y = c(Untreated26, Untreated27, SP_Hp25, SP_Hp26, SP_Hp27), 
                       add.cell.ids = c("Untreated25", "Untreated26", "Untreated27", 
                                        "SP-Hp treated25", "SP-Hp treated26", "SP-Hp treated27"), 
                       project = "Sepsis-Haptoglobin")
pbmc.combined

#Scatter 플랏

plot <- VlnPlot(pbmc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width = 12, repr.plot.height = 6)
plot1 + plot2

pbmc.combined <- subset(pbmc.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
pbmc.combined <- NormalizeData(pbmc.combined, verbose = FALSE)
pbmc.combined <- FindVariableFeatures(pbmc.combined, verbose = FALSE)
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)

saveRDS(pbmc.combined, file = "merging_seurat_data(prototype).rds")

pbmc.combined <- readRDS(file = "d:/scRNAseq/figure/merging_seurat_data(prototype).rds")

write.csv(pbmc.combined@active.ident, file = "d:/cell_population.csv")

plot <- VlnPlot(pbmc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width = 12, repr.plot.height = 6)
plot1 + plot2

prop.table(table(Idents(pbmc.combined)))
table(Idents(pbmc.combined))
View(pbmc.combined)
pbmc.combined$orig.ident <- factor(x = pbmc.combined$orig.ident, levels = c("Untreated", "SP-Hp treated"))
mean(pbmc.combined@meta.data[["nCount_RNA"]])
Idents(pbmc.combined) <- pbmc.combined$orig.ident
plot

plot <- ElbowPlot(pbmc.combined)
options(repr.plot.width = 6, repr.plot.height = 6)
plot

# Clustering--------------------------------------------------------------------
# Clustering--------------------------------------------------------------------
# Clustering--------------------------------------------------------------------

#엘보우 포인트 대략 6-7
pbmc.combined2 <- FindNeighbors(pbmc.combined, dims = 1:8, verbose = FALSE)
pbmc.combined2 <- FindClusters(pbmc.combined2, resolution = 2, verbose = FALSE)
# 일반적으로 resolution 값은 0.1 ~ 1.0 사이의 값을 많이 사용합니다. 
# Resolution 값이 클수록 세분화된 군집을 얻을 수 있기 때문에
# 세포의 종류나 상태 등을 더 세부적으로 파악하고자 할 때는 큰 값이 유용합니다. 
# 반면 작은 값은 대부분의 데이터를 하나의 군집으로 묶어줌으로써 전체적인 데이터 구조를 파악하는 데에 유용할 수 있습니다.
pbmc.combined2 <- RunUMAP(pbmc.combined2, dims = 1:7, verbose = FALSE) # 디멘션 건들면 모양 바뀜. 
# pbmc.combined2 <- RunTSNE(pbmc.combined2, dims = 1:7, verbose = FALSE)

plot <- DimPlot(pbmc.combined2, reduction = "umap", label=TRUE)
options(repr.plot.width = 6, repr.plot.height = 6)
plot

# tplot <- DimPlot(pbmc.combined2, reduction = "tsne", label=TRUE)
# options(repr.plot.width = 6, repr.plot.height = 6)
# tplot
#remove(tplot)

pbmc.combined2$orig.ident <- factor(x = pbmc.combined2$orig.ident, levels = c("Untreated", "SP-Hp treated"))
plot <- DimPlot(pbmc.combined2, reduction = "umap", label=TRUE)
plot

table(pbmc.combined$orig.ident)

pbmc.combined[["cluster"]] <- pbmc.combined@active.ident

pbmc.combined2$cluster

pbmc.combined2 <- readRDS(file = "merging_seurat_data(cluster0_30).rds") #위치 확인
Idents(pbmc.combined2) <- pbmc.combined2$seurat_clusters
pbmc.markers <- FindAllMarkers(pbmc.combined2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top10
DoHeatmap(pbmc.combined2, features = top10$gene) + NoLegend()

anot_genes <- c("PIM2", "RGS1", "CCR7", "FCGR1A", "CD14", "KLRF1", "IL3RA", "LYZ", "CD4", "CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", 
                "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", 
                "CCL2", "S100A9", "HLA-DQA1", "GPR183", "CD103", "CD69", "CD45RA", "CD45RO", "CD62L", "CCR7", "IL7R", "CXCR5", "CD44")

anot_genes2 <- c("LYZ", "MIR155HG", "CCL2", "CST3", 
                 "CD14", "CD4", "FCGR3A", 
                 "CD3D", "SELL", 
                 "CD8A", 
                 "HSPH1", "CD79A", "NKG7",
                 "KLRF1", "GNLY", "CCL5", "GPR183", "CCR7")

anot_genes3 <- c("CD206", "LY6C", "CCL2", 
                 "CD14", "CD4", "FCGR3A", 
                 "CD3D", "SELL", 
                 "CD8A", 
                 "HSPH1", "CD79A",
                 "KLRF1", "GNLY", "CCL5", "GPR183", "CCR7")


View(pbmc.combined)
pbmc.combined@active.ident 

x <- levels(pbmc.combined@active.ident)
levels = c("CD4+ T cells", "CD8+ T cells", "B cells", "Macrophages",
  "Monocytes", "NK cells", "DCs")


pbmc.combined@active.ident <- factor(pbmc.combined@active.ident, levels = rev(c("Macrophages", "Monocytes", "T cells", "T cells", "T cells", "B cells")))
pbmc.combined2 <- pbmc.combined
pbmc.combined2@active.ident[pbmc.combined2@active.ident == "Memory T cells"] <- "CD4+ T cells"
p <- DotPlot(pbmc.combined2, features = anot_genes2, cols = c("lightgrey", "#2ADBA2")) + RotatedAxis()
p + scale_size(range = c(1, 10), breaks = c(0, 20, 40, 60, 80, 100))

pbmc.combined3 <- RenameIdents(
  object = pbmc.combined2,
  `Memory T cells` = "CD4+ T cells",
)

plot <- DimPlot(pbmc.combined, reduction = "umap", label=TRUE)
options(repr.plot.width = 6, repr.plot.height = 6)
plot

pbmc.combined3[["old.ident"]] <- Idents(object = pbmc.combined3)

# 레이블 변경하기
pbmc.combined3 <- RenameIdents(
  object = pbmc.combined3,
  `0` = "CD4+ T cells",
  `1` = "NK cells",
  `2` = "CD4+ T cells",
  `3` = "NK cells",
  `4` = "Monocytes",
  `5` = "CD4+ Memory T cells",
  `6` = "Macrophages",
  `7` = "CD8+ T cells",
  `8` = "B cells", 
  `9` = "B cells", 
  `10` = "DCs", 
  `11` = "CD8+ T cells"
)



pbmc.combined3$seurat_clusters
p <- DimPlot(pbmc.combined3, label = TRUE, group.by = "ident")


pbmc.markers <- FindAllMarkers(pbmc.combined3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers, file = "pbmc_markers.rds")
write.csv(pbmc.markers, file = "pbmc_markers.csv")
pbmc.markers <- readRDS(file = "pbmc_markers.rds")

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10


cells.1 <- WhichCells(pbmc.combined3, idents = 'Monocytes')
cells.2 <- WhichCells(pbmc.combined3, idents = 'Macrophages')
hm1 <- DoHeatmap(pbmc.combined3, cells = cells.1)
hm2 <- DoHeatmap(pbmc.combined3, cells = cells.2)

DoHeatmap(pbmc.combined3, features = top10$gene) + NoLegend()
View(pbmc.combined3)

pbmc.combined3@meta.data$orig.ident

pbmc.combined3[['detail.ident']] <- Idents(object = pbmc.combined2)

#DoMultiBarHeatmap(pbmc.combined3, features = top10$gene, group.by='orig.ident', additional.group.by = 'detail.ident')

plot_heatmap(dataset = pbmc.combined3, 
             markers = pbmc.markers, 
             sort_var = c("seurat_clusters","sample"),
             anno_var = c("seurat_clusters","sample","percent.mt","S.Score","G2M.Score"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("red","orange","yellow","purple","blue","green"), # color vector
                                "Reds",
                                c("blue","white","red"),                            # Three-color gradient
                                "Greens"))


plot_heatmap(dataset = pbmc.combined3, 
             markers = top10,
             sort_var = c("seurat_clusters","sample"),
             anno_var = c("seurat_clusters","sample","percent.mt","S.Score","G2M.Score"),
             anno_colors = list("Set2",                                             # RColorBrewer palette
                                c("red","orange","yellow","purple","blue","green"), # color vector
                                "Reds",
                                c("blue","white","red"),                            # Three-color gradient
                                "Greens"))

# Cell type identification------------------------------------------------------
# Cell type identification------------------------------------------------------
# Cell type identification------------------------------------------------------

# 각각의 클러스터들이 어떤 세포인지 알아내는 작업을 합니다.Seurat의 FindAllMarkers()함수를 사용하면 나머지 모든 세포와 비교해 클러스터에 대한 유전자 마커를 찾을 수 있습니다.
# Marker gene 찾기

#DEG 확인용
#markers <- FindAllMarkers(pbmc.combined2, only.pos = FALSE, verbose = FALSE, logfc.threshold = 0.01)

markers <- FindAllMarkers(pbmc.combined2, only.pos = T, verbose = FALSE)
saveRDS(markers, file = "markers.rds")

write.csv(markers, file = "Allset_marker(resolution2).csv") # 결과를 csv 파일로 저장
markers %>% head()

saveRDS(pbmc.combined2, file = "merging_seurat_data(cluster0_30).rds")
setwd("C:/Users/User/Documents/")
pbmc.combined2 <- readRDS(file = "merging_seurat_data(cluster0_30).rds") #위치 확인
View(pbmc.combined2)
pbmc.combined2$seurat_clusters
p <- DimPlot(pbmc.combined2, label = TRUE, group.by = "ident")

remove(markers)
# Marker gene heatmap-----------------------------------------------------------
# Marker gene heatmap-----------------------------------------------------------
# Marker gene heatmap-----------------------------------------------------------

print(pbmc.combined2[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
pbmc.combined2 <- ScaleData(pbmc.combined2, verbose = FALSE)


p <- DoHeatmap(object = pbmc.combined2, ) + NoLegend()

options(repr.plot.width = 7, repr.plot.height = 7)
p

# ref1 <- HumanPrimaryCellAtlasData()
# ref2 <- NovershternHematopoieticData()
# ref3 <- DatabaseImmuneCellExpressionData()
# ref4 <- MonacoImmuneData()

DefaultAssay(pbmc.combined2) <- "RNA"
combine.all.sce <- as.SingleCellExperiment(pbmc.combined2)
pred.Result <- SingleR(test = combine.all.sce, ref = ref1, assay.type.test=1, labels = ref1$label.fine)
tab <- table(cluster = pbmc.combined2$seurat_clusters, label = pred.Result$labels)
pheatmap::pheatmap(log10(tab+10), cluster_cols = F)


# Labeling ---------------------------------------------------------------------
# Labeling ---------------------------------------------------------------------
# Labeling ---------------------------------------------------------------------
# 참조를 위해 이전 ID 클래스(클러스터 레이블)를 저장합니다.
pbmc.combined2[["old.ident"]] <- Idents(object = pbmc.combined2)
Idents(object = pbmc.combined2) <- pbmc.combined2$seurat_clusters

# 레이블 변경하기
pbmc.combined2 <- RenameIdents(
  object = pbmc.combined2,
  `0` = "CD4+ T cells",
  `1` = "NK cells",
  `2` = "CD4+ T cells",
  `3` = "NK cells",
  `4` = "Monocytes",
  `5` = "CD4+ T cells",
  `6` = "Macrophages",
  `7` = "CD8+ T cells",
  `8` = "B cells", 
  `9` = "B cells", 
  `10` = "DCs", 
  `11` = "CD8+ T cells"
)

pbmc.markers2 <- FindAllMarkers(pbmc.combined2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Idents(pbmc.combined2) <- factor(x = Idents(pbmc.combined2), levels = c("Untreated", "SP-Hp treated"))
pbmc.combined$orig.ident <- factor(x = pbmc.combined$orig.ident, levels = c("Untreated", "SP-Hp treated"))
levels = c("DCs", "NK cells", "CD4+ T cells", "CD8+ T cells", "B cells", "Monocytes", "Macrophages")
Idents(pbmc.combined2) <- factor(x = Idents(pbmc.combined2), levels = levels)
DotPlot(pbmc.combined2, features = anot_genes2, cols = c("lightgrey", "#E7211A"), dot.scale = 10) + RotatedAxis()

p <- DimPlot(pbmc.combined2, label = T, reduction = "tsne", split.by = "orig.ident")
options(repr.plot.width = 7, repr.plot.height = 7)
p

saveRDS(pbmc.combined2, file = "Preprocessing_scRNAseq_data.rds")
saveRDS(pbmc.combined2, file = "Allset_seurat0.5(HumanPrimaryCellAtlasData).rds")



###################### Save point ##############################################
pbmc.combined2 <- readRDS(file = "d:/scRNAseq/figure/Allset_seurat0.5(HumanPrimaryCellAtlasData).rds")
pbmc.combined2$orig.ident <- factor(x = pbmc.combined2$orig.ident, levels = c("Untreated", "SP-Hp treated"))

levels(Idents(pbmc.combined2))

pbmc.combined2 <- RenameIdents(
  object = pbmc.combined2,
  `CD4+ T cells` = "CD4+ T cells",
  `NK cells` = "NK cells",
  `Monocytes` = "Monocytes",
  `CD4+ Memory T cells` = "CD4+ T cells",
  `Macrophages` = "Macrophages",
  `CD8+ T cells` = "CD8+ T cells",
  `Naive B cells` = "B cells", 
  `Memory B cells` = "B cells", 
  `DCs` = "DCs", 
)

p <- DimPlot(pbmc.combined2, label = TRUE, group.by = "ident")


eset <- cbind(pbmc.combined2@meta.data, pbmc.combined2@active.ident)
write.csv(eset, file = "scRNAseq_cell_proportion.csv")

eset2 <- pbmc.combined2@assays[["RNA"]]@scale.data
saveRDS(eset2, file = "scRNAseq_exp.rds")

p1 <- VlnPlot(UnMono_SSMacro, 
        features = c("NLRP3", "IL18", "CASP3", "CASP1", "IL6", "TNF", "CCL4", "CCL3", "CCL2", "CXCL8", "CLEC4E", "CD38"), 
        group.by = "orig.ident", 
        ncol = 4, 
        #pt.size = 0, 
        combine = FALSE, 
        y.max = 8)

wrap_plots(plots = p1, ncol = 4)


head(eset)


#BiocManager::install('dittoSeq')

library(dittoSeq)


RidgePlot(pbmc.combined2, features = 'IL1B', split.by = "orig.ident")
VlnPlot(pbmc.combined2, "CD19")

dittoBarPlot(
  object = pbmc.combined2,
  var = "RNA_snn_res.0.5",
  group.by = "orig.ident") +
  scale_x_discrete(limit=c("Untreated", "SP-Hp treated"))

prop.table(table(Idents(pbmc.combined2)))
table(Idents(pbmc.combined2))

VlnPlot(UnMono_SSMacro, 
        features = c("TNF"), 
        group.by = "orig.ident", 
        ncol = 4, 
        y.max = 8
        )


VlnPlot(UnMono_SSMacro, 
        features = c("CX3CR1", "OR6B2"), 
        group.by = "orig.ident", 
        ncol = 4, 
)


# Gene expression---------------------------------------------------------------
# Gene expression---------------------------------------------------------------
# Gene expression---------------------------------------------------------------

FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("IL1B"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("IL6"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("TNF"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("CXCL10"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("CCL4"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("CCL3"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("CXCL8"), split.by = "orig.ident")
FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("CCL2"), split.by = "orig.ident")
VlnPlot(pbmc.combined2, features = c("IL1B", "IL6", "TNF", "CD38", "CCL4", "CCL3", "CXCL8", "CCL2"), split.by = "orig.ident", split.plot = T)

FeaturePlot(pbmc.combined2,  label = T, reduction = "umap", features = c("CCL2"), split.by = "groups")
pbmc.combined2$orig.ident
pbmc.combined2$old.ident

View(pbmc.combined2)


# 조건에 따른 차등발현 유전자 확인하기
idx <- which(pbmc.combined2@meta.data$orig.ident == "Untreated")
UnMM <- pbmc.combined2[,idx]
UnMM <- subset(UnMM, ident = c("Monocytes"))
Idents(UnMM) <- "Un"
# avg_UnMM <- log1p(AverageExpression(UnMM, verbose = FALSE)$RNA)
# avg_UnMM$gene <- rownames(avg_UnMM)

idx <- which(pbmc.combined2@meta.data$orig.ident == "SP-Hp treated")
SSMM <- pbmc.combined2[,idx]
SSMM <- subset(SSMM, ident = c("Macrophages"))
Idents(SSMM) <- "SP"
# avg_SSMM <- log1p(AverageExpression(SSMM, verbose = FALSE)$RNA)
# avg_SSMM$gene <- rownames(avg_SSMM)

df= pbmc[["SCT"]]@scale.data



UnMono_SSMacro <- merge(UnMM, y = SSMM, 
                       add.cell.ids = c("Un-Monocytes", "SP-Hp-Macrophages"), 
                       project = "Sepsis-Haptoglobin")

saveRDS(UnMono_SSMacro, file = "UnMONO_vs_SPMacro.rds")


#################### Save point ################################################
UnMono_SSMacro <- readRDS(file = "d:/scRNAseq/figure/UnMONO_vs_SPMacro.rds")

VlnPlot(UnMono_SSMacro, features = c("IL1B", "IL6", "TNF", "CCL4", "CCL3", "CCL2", "CXCL8", "CLEC4E", "CX3CL1"), group.by = "orig.ident", ncol = 4)

corr <- FindVariableFeatures(SSMM, selection.method = "vst", nfeatures = 2000)
View(corr)
sd <- corr@assays$RNA@scale.data
View(sd)
length(colnames(sd))
write.csv(sd, file = "SSMM_scale_data.csv")


###### Correlation -------------------------------------------------------------
IF <- read.csv(file = "d:/IF_response.csv")
IF$Symbol <- toupper(IF$Symbol)
IFgenes <- IF$Symbol
tsd2 <- tsd[, colnames(tsd) %in% IFgenes]
dim(tsd)
dim(tsd2)
colnames(tsd2)
tsd <- as.data.frame(t(sd))
cor.test(tsd$CLEC4E, tsd$IL1B)

target <- "CLEC4E"
target3 <- tsd[, "CLEC4E"]




test <- cor.test(target3, tsd2[, 1])
View(test)
test$p.value
test$estimate
count <- dim(tsd2)[2]
EP <- cbind(Rho = test$estimate, `p value` = test$p.value)
colnames(tsd2)[1]

tsd3 <- cbind(target3, tsd2)
colnames(tsd3)[1] <- "CLEC4E"

for(i in 1:count) {
  temp <- cor.test(target3, tsd2[, i])
  EP <- cbind(Rho = temp$estimate, `p value` = temp$p.value)
  rownames(EP) <- colnames(tsd2)[i]
  
  if(i == 1){
    data <- EP 
  } else {
    data <- rbind(data, EP)
  }
}

data <- as.data.frame(data)
str(data)

data2 <- data  %>%
  filter(!is.na(Rho))

if(round(data2["NLRP3",2], 4) == 0) {
  pval <- "0.000"
} else {
  pval <- round(data2["NLRP3",2], 4)
}

rlabel <- paste0("Rho = ", signif(round(data2["NLRP3",1], 4)), '\n', "p = ", pval, " ")



ggplot(tsd3, aes(x = CLEC4E, y = IL1B)) +
  geom_point(col = "black", size=2, alpha = 0.2) +
  theme_bw() +
  xlab("CLEC4E") +
  ylab("IL1B") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  
#   annotate("text", x=mean(range(faithful$eruptions)), y=-Inf, vjust=-0.4, label= "Bottom middle") 

ggplot(tsd3, aes(x = CLEC4E, y = CASP1)) +
  geom_point(col = "black", size=2, alpha = 0.2) +
  theme_bw() +
  xlab("CLEC4E") +
  ylab("CASP1") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  
   annotate("text", x=mean(range(faithful$eruptions)), y=-Inf, vjust=-0.4, label= "Bottom middle") 

ggplot(tsd3, aes(x = CLEC4E, y = NLRP3)) +
  geom_point(col = "black", size=2, alpha = 0.2) +
  theme_bw() +
  xlab("CLEC4E") +
  ylab("NLRP3") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  
#   annotate("text", x=mean(range(faithful$eruptions)), y=-Inf, vjust=-0.4, label= "Bottom middle") 

data <- cbind(Symbol = gene$sym, data)

write.csv(data, file = "CLEC4E_corr.csv")
plot(target3, tsd2$IL1B)




############ Heatmap -----------------------------------------------------------
library(pheatmap)

head(tsd3)
data3 <- data2[data2$Rho > 0 & data2$`p value` < 0.05, ]
head(data2)
length(rownames(data3))

cor_tsd3 <- cor(tsd3[, c("CLEC4E", rownames(data3))])
dim(cor_tsd3)



colfunc <- colorRampPalette(c("blue", "white"))
colfunc2 <- colorRampPalette(c("white", "red"))

postscript(file="Mmu_IF_Heatmap.eps", onefile=FALSE, horizontal=FALSE, width = 1000, height = 3000)

x <- pheatmap(cor_tsd3, 
              color = c(colfunc(100), colfunc2(100)),
              border_color = "white",
              #cutree_cols = 3,
              #cutree_rows = 1,
              cellwidth = 20,
              cellheight = 20,
              cluster_rows = F,
              cluster_cols = F,
              show_colnames = T, 
              fontsize = 8,
              #main = "Inflammatory Cytokine/Chemokine Genes Heatmap \n (Z-score)",
              #scale = "row",
              breaks = seq(from = -1, to = 1, length.out = 200),
              angle_col = "0")
#  width = 800,
#  height = 600,
#gaps_col = c(3,6),
#gaps_row = tempb,
#legend_breaks =  c(-1, 0, 1),
#legend_labels = -2:2,

dev.off()
############ Heatmap -----------------------------------------------------------





length(UnMM$old.ident)
length(SSMM$old.ident)

mf <- corr@assays$RNA@meta.features
View(mf)
length(colnames(mf))

ct <- corr@assays$RNA@counts
View(ct)
colnames(ct)


deg_UnMono_SSMacro <- FindMarkers(UnMono_SSMacro, ident.1 = "SP", ident.2 = "Un", verbose = FALSE)
write.csv(deg_UnMono_SSMacro, file = "deg_UnMono_SSMacro.csv")
dim(deg_UnMono_SSMacro)

genes.to.label = c("CX3CL1", "IL1B", "IL6", "TNF", "CCL4", "CCL3", "CCL2", "CXCL10", "CXCL8")
p1 <- ggplot(avg_UnMM, aes(Untreated, `SP-Hp treated`)) + geom_point() + ggtitle("Untreated_Monocytes & Macorphages")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)

pbmc.combined2$celltype.stim <- paste(Idents(pbmc.combined2), immune.combined$stim, sep = "_")

levels(pbmc.combined2)
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")


DimPlot(pbmc.combined2, reduction = "umap", group.by = "orig.ident")

markers.to.plot <- c("IL1B", "IL6", "TNF", "CCL4", "CCL3", "CCL2", "CXCL10", "CX3CR1")
DotPlot(pbmc.combined2, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "orig.ident") + RotatedAxis()


saveRDS(pbmc.combined2)

























VlnPlot(seurat_obj, features = "IL6")
#VlnPlot(seurat_obj, features = "IL1B", idents = "Monocytes")


saveRDS(seurat_obj, file = "HC25.rds")
View(seurat_obj$nFeature_RNA)

seurat_obj <- readRDS(file = "D:/scRNAseq/HC25.rds")
seurat_obj2 <- readRDS(file = "D:/scRNAseq/SS25.rds")

FeaturePlot(seurat_obj2,  reduction = "umap", features = c("IL6"))

VlnPlot()

VlnPlot(temp, features = "IL1B", group.by = "orig.ident", idents = "Monocytes", split.by = "orig.ident")
VlnPlot(seurat_obj, features = "IL6") + VlnPlot(seurat_obj2, features = "IL6")
VlnPlot(seurat_obj, features = "TNF") + VlnPlot(seurat_obj2, features = "TNF")
#VlnPlot(seurat_obj, features = "IL10") + VlnPlot(seurat_obj2, features = "IL10")
#VlnPlot(seurat_obj, features = "IL18") + VlnPlot(seurat_obj2, features = "IL18")
VlnPlot(seurat_obj, features = "CXCL10") + VlnPlot(seurat_obj2, features = "CXCL10")
VlnPlot(seurat_obj, features = "CCL4") + VlnPlot(seurat_obj2, features = "CCL4")
VlnPlot(seurat_obj, features = "CCL3") + VlnPlot(seurat_obj2, features = "CCL3")
VlnPlot(seurat_obj, features = "CXCL8") + VlnPlot(seurat_obj2, features = "CXCL8")
VlnPlot(seurat_obj, features = "CCL2") + VlnPlot(seurat_obj2, features = "CCL2")
#VlnPlot(seurat_obj, features = "IFNG") + VlnPlot(seurat_obj2, features = "IFNG")
#VlnPlot(seurat_obj, features = "CSF3") + VlnPlot(seurat_obj2, features = "CSF3")


###### Merge -------------------------------------------------------------------

pbmc.combined <- merge(seurat_obj, y = seurat_obj2, add.cell.ids = c("Untreated", "SP-Hp"), project = "Sepsis")
pbmc.combined

table(pbmc.combined$orig.ident)

pbmc.combined <- RunUMAP(object = pbmc.combined, dims=1:6)
pbmc.combined <- RunTSNE(object = pbmc.combined, dims=1:6)
FeaturePlot(pbmc.combined,  reduction = "tsne", features = c("IL6"))

DefaultAssay(pbmc.combined) <- "RNA"  # default assay is RNA
pbmc.combined <- NormalizeData(pbmc.combined, verbose = FALSE)
pbmc.combined <- FindVariableFeatures(pbmc.combined, verbose = FALSE)
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:10, verbose = FALSE)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:10, verbose = FALSE)

p <- DimPlot(pbmc.combined, label = TRUE, group.by = "orig.ident")
options(repr.plot.width = 7, repr.plot.height = 7)
p

DefaultAssay(pbmc.combined) <- "RNA"
combine.all.sce <- as.SingleCellExperiment(pbmc.combined)
pred.Result <- SingleR(test = combine.all.sce, ref = ref1, assay.type.test=1, labels = ref1$label.fine)

DimPlot(pbmc.combined, reduction = "umap", label = T, repel = T)
tab <- table(cluster = pbmc.combined$seurat_clusters, label = pred.Result$labels)
pheatmap::pheatmap(log10(tab+10), cluster_cols = F)

pbmc.combined[["old.ident"]] <- Idents(object = pbmc.combined)

# 레이블 변경하기
pbmc.combined <- RenameIdents(
  object = pbmc.combined,
  `0` = "CD4+ T cells",
  `1` = "NK cells",
  `2` = "CD4+ T cells",
  `3` = "CD8+ T cells",
  `4` = "CD4+ T cells",
  `5` = "Naive B cells",
  `6` = "Monocytes",
  `7` = "Memory B cells",
  `8` = "Monocytes", 
  `9` = "CD8+ T cells", 
  `10` = "DCs"
)

p <- DimPlot(pbmc.combined, label = TRUE, reduction = "tsne", group.by = "orig.ident")

options(repr.plot.width = 7, repr.plot.height = 7)
p
