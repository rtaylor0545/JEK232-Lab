eset <- read.csv(file = "d:/Newseq.csv")
eset[1:5, 1:5]
IF <- read.csv(file = "d:/GO_term/Hsa_Inflammation_genes.csv")
IF2 <- IF$HM_Hsa_Inflammatory.response[IF$HM_Hsa_Inflammatory.response != ""]
eset_IF <- eset[eset$Gene_Symbol %in% IF2, ]
dim(eset_IF)
eset_IF2 <- eset_IF[,c(1,11:dim(eset_IF)[2])]
rownames(eset_IF2) <- eset_IF2$Gene_Symbol
eset_IF3 <- eset_IF2[,-1]

eset_IF2order <- eset_IF[order(-eset_IF$G_O.fc),]
eset_IF2order$Gene_Symbol

eset_IF2order2 <- eset_IF2order[eset_IF2order$G_B.fc > 1, ]
eset_IF2order3 <- eset_IF2order2[eset_IF2order2$G_O.fc > 1, ]
dim(eset_IF2order3)
eset_IF2order3$Gene_Symbol

write.csv(eset_IF2order3, file = "IF_DEG.csv")

eset_IF4 <- eset_IF3[eset_IF2order3$Gene_Symbol,]

library("pheatmap")
library("RColorBrewer")

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Blues")[7]))(19))

a <- pheatmap(eset_IF4, 
              scale = "row", 
              legend = T, 
              color = cols,
              border_color = "black", 
              gaps_col = c(12, 31),
              clustering_distance_rows = "euclidean",
              fontsize = 5,
              #gaps_row = 1:dim(mf7)[1],
              cellwidth = 6, 
              cellheight = 6,
              #legend_breaks = c(0.2, 0.4, 0.6),
              breaks = seq(-1.5, 1.5, 3/19),
              show_colnames = F,
              show_rownames = T,
              cluster_rows = F,
              cluster_cols = F)
