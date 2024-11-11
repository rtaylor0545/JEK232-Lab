setwd("C:/Users/User/Documents")
#BiocManager::install('Seurat')
#BiocManager::install('Matrix', version = '1.6-4')

library(Seurat)

pbmc.combined2 <- readRDS(file = "merging_seurat_data(cluster0_30).rds")
plot <- DimPlot(pbmc.combined2, reduction = "umap", label=TRUE)

Idents(pbmc.combined2) <- "orig.ident"

DimPlot(pbmc.combined2, reduction = "umap", label=F, group.by = "orig.ident")
DimPlot(pbmc.combined2, reduction = "umap", label=F, split.by = "orig.ident")

# 특정 클러스터만 추리기
Monocytes <- subset(pbmc.combined2, ident = c("27", "2", "19"))
Macrophages <- subset(pbmc.combined2, ident = c("29", "18", "7"))

test <- subset(pbmc.combined2, ident = c("27", "2", "19", "29", "18", "7"))

test$celltype <- "Monocyte"
test$celltype[test$seurat_clusters %in% c("29", "18", "7")] <- "Macrophage"
Idents(test) <- "celltype"

test@active.ident <- 
  factor(test@active.ident, 
         levels = c("7", "18", "29", "19", "2", "27"))

VlnPlot(test, features = c("FUT3")) + NoLegend() + coord_flip() 
ggsave(filename = paste0("FUT3", ".pdf"), dpi = 72, width = 150, height = 350, units = 'px') #12)

VlnPlot(test, features = c("FUT4")) + NoLegend() + coord_flip() 
ggsave(filename = paste0("FUT4", ".pdf"), dpi = 72, width = 150, height = 350, units = 'px') #12)

VlnPlot(test, features = c("FUT10")) + NoLegend()+ coord_flip() 
ggsave(filename = paste0("FUT10", ".pdf"), dpi = 72, width = 150, height = 350, units = 'px') #12)

VlnPlot(test, features = c("FUT11")) + NoLegend()+ coord_flip() 
ggsave(filename = paste0("FUT11", ".pdf"), dpi = 72, width = 150, height = 350, units = 'px') #12)

write.table(test@active.ident, file='Convert_UMI_Label.tsv', quote=FALSE, sep='\t', col.names = TRUE)
write.table(test@assays[["RNA"]]@counts, file='Gene_Count_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)

getwd()

test.markers <- FindAllMarkers(test, only.pos = TRUE)
test.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

test.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

write.csv(test.markers, file = "testmarker.csv")

plot <- DimPlot(test, reduction = "umap", label=TRUE)

DoHeatmap(test, features = top10$gene) + NoLegend()

mono.markers <- FindAllMarkers(Monocytes, only.pos = T)
macro.markers <- FindAllMarkers(Macrophages, only.pos = T)

library(dplyr)

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

macro.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> macro.top10

DoHeatmap(Macrophages, features = macro.top10$gene) + NoLegend()
p <- DotPlot(Macrophages, features = macro.top10$gene, cols = c("lightgrey", "#2ADBA2")) + RotatedAxis()

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

