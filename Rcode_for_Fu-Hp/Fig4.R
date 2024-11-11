### Fig4d ###

setwd("C:/Users/User/Documents")
#BiocManager::install('Seurat')
#BiocManager::install('Matrix', version = '1.6-4')

library(Seurat)

pbmc.combined <- readRDS(file = "D:/scRNAseq/figure/Allset_seurat0.5(HumanPrimaryCellAtlasData).rds")
plot <- DimPlot(pbmc.combined, reduction = "umap", label=TRUE)

pbmc.combined <- RenameIdents(
  object = pbmc.combined,
  `CD4+ T cells` = "CD4+ T cells",
  `CD4+ Memory T cells` = "CD4+ T cells",
  `CD8+ T cells` = "CD8+ T cells",
  `NK cells` = "NK cells",
  `Monocytes` = "Monocytes",
  `Macrophages` = "Macrophages",
  `Naive B cells` = "B cells",
  `Memory B cells` = "B cells",
  `DCs` = "DCs"
)

pbmc.combined2$orig.ident

pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE)
View(pbmc.markers)
write.csv(pbmc.markers, "d:/scRNAseq(Fuc-Hp).csv")

DotPlot(pbmc.combined, features = c("CD80", "CD86", "FCGR1A", "HLA-A", "CD200R1", "MRC1", "CD163"), scale.max = 40)
saveRDS(pbmc.markers, file = "d:/scRNAseq(Fuc-Hp_markers).rds")
pbmc.markers <- readRDS(file = "d:/scRNAseq(Fuc-Hp_markers).rds")


library(stringr)

head(names(pbmc.combined2$orig.ident))
tail(names(pbmc.combined2$orig.ident))

pbmc.combined2$celltype <- "Untreated(1)"
pbmc.combined2$celltype[(str_detect(names(pbmc.combined2$orig.ident), "Untreated26"))] <- "Untreated(2)"
pbmc.combined2$celltype[(str_detect(names(pbmc.combined2$orig.ident), "Untreated27"))] <- "Untreated(3)"

pbmc.combined2$celltype[(str_detect(names(pbmc.combined2$orig.ident), "SP-Hp treated25"))] <- "SS-Hp(1)"
pbmc.combined2$celltype[(str_detect(names(pbmc.combined2$orig.ident), "SP-Hp treated26"))] <- "SS-Hp(2)"
pbmc.combined2$celltype[(str_detect(names(pbmc.combined2$orig.ident), "SP-Hp treated27"))] <- "SS-Hp(3)"

pbmc.combined2$Cell_annot <- Idents(pbmc.combined2)

plot <- DimPlot(pbmc.combined2, reduction = "umap", label=TRUE)

library(ggplot2)
library(RColorBrewer)

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Blues")[7]))(7))

levels <- levels(factor(pbmc.combined2@meta.data$celltype))
le1 <- rev(levels[1:3])
le2 <- rev(levels[4:6])
pbmc.combined2@meta.data$celltype <- factor(pbmc.combined2@meta.data$celltype, levels= c(le1, le2))

ggplot(pbmc.combined2@meta.data, aes(x=celltype, fill=Cell_annot)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols) + coord_flip()

Idents(pbmc.combined2) <- pbmc.combined2$orig.ident
pbmc.combined2@meta.data$orig.ident <- factor(pbmc.combined2@meta.data$orig.ident, levels= c("Untreated", "SP-Hp treated"))
plot <- DimPlot(pbmc.combined2, reduction = "umap", label=F, split.by = "orig.ident") + NoLegend()

Idents(pbmc.combined2) <- pbmc.combined2$Cell_annot
pbmc.combined3 <- subset(pbmc.combined2, ident = c("Monocytes", "Macrophages"))

cols <-  (colorRampPalette(c( 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Reds")[7]))(5))


pbmc.combined3@meta.data$orig.ident <- factor(pbmc.combined3@meta.data$orig.ident, levels= c("SP-Hp treated", "Untreated"))
ggplot(pbmc.combined3@meta.data, aes(x=orig.ident, fill=Cell_annot)) + geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(values=cols)

pbmc.combined2 <- readRDS(file = "merging_seurat_data(cluster0_30).rds")
pbmc.combined3 <- subset(pbmc.combined2, ident = c("27", "2", "19", "29", "18", "7"))
pbmc.combined3@active.ident <- factor(x = pbmc.combined3@active.ident, levels = rev(c("29", "18", "7", "27", "19", "2")))

DotPlot(pbmc.combined, features = c("CD14","CD64", "CD40","CD32", 
                                    "CD163", "CD206", "CD80", "CD86", 
                                    "FCGR1A", "HLA-A", "CD200R1", "MRC1", "CD68", "ADGRE1"), scale.max = 50)

DotPlot(pbmc.combined, features = c("PTX3","CCL15", "IL12B","CD14", 
                                    "CD68", "CD163"), scale.max = 25)

DotPlot(pbmc.combined, features = c("IL6", "CCL15", "IL12B", "TNF", "IL15RA", "BCL2A1"), scale.max = 10)
VlnPlot(pbmc.combined, features = c("PTX3", "CXCL5", "CCL15", "CCL23", "CXCL6", "IL20", "IL12B", "CXCL1", "TNFAIP6"))
DotPlot(pbmc.combined, features = c("PTX3", "CXCL5", "CCL15", "CCL23", "CXCL6", "IL20", "IL12B", "CXCL1", "TNFAIP6"), scale.max = 25)
DotPlot(pbmc.combined, features = c("PTX3", "CCL15","IL12B", "IL1B", "TNF", "IL6", "CCL20", 
                                    "SLC21A15", "INDO", "HSD11B1", "SPHK1", "PFKFB3", "PSMA2"), scale.max = 25)

VlnPlot(pbmc.combined3, features = c("CXCL1", "CXCL8", "CXCL5", "SERPINB2", "IL1B", "IL23A"))
VlnPlot(pbmc.combined3, features = c("CXCL6", "CXCL5", "F3", "SERPINB2", "IL1B", "CXCL3"))

library(dplyr)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) -> top10
DoHeatmap(pbmc.combined, features = top10$gene) + NoLegend()


DotPlot(pbmc.combined3, features = c("IL1B", "CCL3", "TNF", "CXCL16", "CXCL6", "CXCL5"), scale.max = 25)

pbmc.combined4 <- subset(pbmc.combined, ident = c("Monocytes", "Macrophages"))

plot <- DimPlot(pbmc.combined3, reduction = "umap", label=T) + NoLegend()

Idents(pbmc.combined4)

deg_pbmc.combined4 <- FindMarkers(pbmc.combined4, ident.1 = "Monocytes", ident.2 = "Macrophages", verbose = FALSE)
write.csv(deg_pbmc.combined4, file = "deg_pbmc.combined4.csv")
dim(deg_UnMono_SSMacro)


# 특정 클러스터만 추리기
Monocytes <- subset(pbmc.combined2, ident = c("27", "2", "19"))
Macrophages <- subset(pbmc.combined2, ident = c("18", "7", "29"))

test <- subset(pbmc.combined2, ident = c("27", "2", "19", "29", "18", "7"))


test@meta.data$seurat_clusters <- factor(test@meta.data$seurat_clusters, levels= c("2", "19", "27", "7", "18", "29"))
Idents(test) <- test$seurat_clusters

RidgePlot(test, features = c("CXCL2", "CXCL3", "TNF", "CCL3", "CXCL8", "CCL4", "IL1B", "IL1A", "IL6", "CCL20", "CXCL1"), ncol = 1) + NoLegend()
RidgePlot(test, features = c("CLEC4E"), ncol = 1) + NoLegend()

test.markers <- FindAllMarkers(Macrophages, only.pos = TRUE)
test.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

test.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

write.csv(test.markers, file = "MAC_testmarker.csv")

Macrophages@active.ident <- factor(x = Macrophages@active.ident, levels = rev(c("18", "7", "29")))

MAC_genes <- c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "MMP9", 
               "IL1B", "CLEC4E", "CCL1", "CCL2", "CCL3", "CCL7", 
               "ALOX5AP", "MMP14", "NFKBIZ", "PKM", "TNF")
p <- DotPlot(Macrophages, features = MAC_genes, cols = c("lightgrey", "#2ADBA2"), dot.scale = 15) + RotatedAxis()
p 

DotPlot(Macrophages, features = c("MMP8", "CLEC4E"), cols = c("lightgrey", "#2ADBA2"), dot.scale = 15) + RotatedAxis()
p + scale_size(range = c(1, 10), breaks = c(0, 25, 50, 100))
RidgePlot(Macrophages, features = "MMP8", ncol = 2)




p <- DotPlot(Macrophages, features = MAC_genes, cols = c("lightgrey", "#2ADBA2")) + RotatedAxis()
p + scale_size(range = c(1, 10), breaks = c(0, 25, 50))


plot <- dotplot(test, reduction = "umap", label=TRUE)

DoHeatmap(test, features = top10$gene) + NoLegend()

mono.markers <- FindAllMarkers(Monocytes, only.pos = T)
macro.markers <- FindAllMarkers(Macrophages, only.pos = T)

mono.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> mono.top10

macro.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 100) %>%
  ungroup() -> macro.top100

plot <- DimPlot(Macrophages, reduction = "umap", label=TRUE)

DoHeatmap(Macrophages, features = macro.top100$gene) + NoLegend()

vplot1 <- VlnPlot(Macrophages, features = c("IL1A", "IL1B", 
                                            "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL10", 
                                            "CCL1", "CCL20", "CCL23", 
                                            "CLEC4E", "CLEC5A"), ncol = 7)

vplot2 <- VlnPlot(Monocytes, features = c("IL1A", "IL1B", 
                                          "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL10", 
                                          "CCL1", "CCL20", "CCL23", 
                                          "CLEC4E", "CLEC5A"), ncol = 7)


test@active.ident <- factor(x = test@active.ident, levels = c("27", "2", "19", "29", "18", "7"))


vplot1 <- VlnPlot(test, features = c("CXCL5"), pt.size = 0, ncol = 1)
vplot2 <- VlnPlot(test, features = c("CXCL1"), pt.size = 0, ncol = 1)
vplot3 <- VlnPlot(test, features = c("CXCL3"), pt.size = 0, ncol = 1)
vplot4 <- VlnPlot(test, features = c("IL1B"), pt.size = 0, ncol = 1)
vplot5 <- VlnPlot(test, features = c("CCL3"), pt.size = 0, ncol = 1)
vplot6 <- VlnPlot(test, features = c("CLEC4E"), pt.size = 0, ncol = 1)

