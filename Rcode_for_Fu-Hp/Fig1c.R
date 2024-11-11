### Fig.1c Heatmap ###

genes <- read.csv(file = "d:/GO_term/IF_TNF.csv")
NCC <- read.csv(file = "d:/NewCNU_Cohort.csv", row.names = 1) 

SS <- NCC[,1:24]
HC <- NCC[,25:36]

SS2 <- SS[rownames(SS) %in% c(genes$IF_Gene, genes$TNFA_Gene),]
HC2 <- HC[rownames(HC) %in% c(genes$IF_Gene, genes$TNFA_Gene),]

dim(SS2)
dim(HC2)

NCC2 <- cbind(HC2, SS2)

overlapped <- intersect(genes$IF_Gene, genes$TNFA_Gene)
IFonly <- genes$IF_Gene[!(genes$IF_Gene %in% overlapped)]
TNFAonly <- genes$TNFA_Gene[!(genes$TNFA_Gene %in% overlapped)]

NCC3 <- rbind(NCC2[IFonly,], NCC2[overlapped,], NCC2[TNFAonly,])

FC <- read.csv(file = "d:/SSvsHC.csv", row.names = 1)
FC_IF <- FC[IFonly,]
FC_OV <- FC[overlapped,]
FC_TNF <- FC[TNFAonly,]

FC_IF2 <- FC_IF[(FC_IF$logFC > 0 & FC_IF$PValue < 0.05), ]
FC_OV2 <- FC_OV[(FC_OV$logFC > 0 & FC_OV$PValue < 0.05), ]
FC_TNF2 <- FC_TNF[(FC_TNF$logFC > 0 & FC_TNF$PValue < 0.05), ]

FC_IF3 <- FC_IF2[order(FC_IF2$PValue),]
FC_IF4 <- FC_IF3[1:20,]

FC_OV3 <- FC_OV2[order(FC_OV2$PValue),]
FC_OV4 <- FC_OV3[1:10,]

FC_TNF3 <- FC_TNF2[order(FC_TNF2$PValue),]
FC_TNF4 <- FC_TNF3[1:20,]

final_genes <- c(rownames(FC_IF4), rownames(FC_OV4), rownames(FC_TNF4))
NCC3 <- NCC2[final_genes,]

NCC3[1:3, 1:3]

library("pheatmap")
library("RColorBrewer")

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Blues")[7]))(19))

dev.off()

length(IFonly)
length(overlapped)
length(TNFAonly)

a <- pheatmap(NCC3, 
              scale = "row", 
              legend = T, 
              color = cols,
              border_color = "black", 
              #gaps_col = c(12, 31),
              clustering_distance_rows = "euclidean",
              fontsize = 8,
              #gaps_row = 1:dim(mf7)[1],
              cellwidth = 5, 
              cellheight = 10,
              legend_breaks = c(-1, 0, 1),
              breaks = seq(-1, 1, 2/19),
              show_colnames = F,
              show_rownames = T,
              cluster_rows = F,
              cluster_cols = F, 
              angle_col = 0, 
              family = "Arial Narrow")
