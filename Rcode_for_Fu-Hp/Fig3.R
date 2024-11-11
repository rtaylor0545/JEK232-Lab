library("pheatmap")
library("RColorBrewer")
library(ggplot2)
library('ggrepel')
library(edgeR)

### edgeR로 DEG 구하기 ###
{
  NCC <- read.csv(file = "c:/Users/User/Desktop/01_Sepsis_Hp_Project/Fig.3/HP-treated_PBMC_RNAseq.csv")  ### Kallisto data (이미 노말라이제이션까지 다 진행됨)
  NCC[1:5, 1:5]
  rownames(NCC) <- NCC$Gene
  NCC <- NCC[,-1]
  NCC[1:5, 1:4]
  
  group <- factor(c(rep('HC', 2), rep('SS', 2)))
  y <- DGEList(counts = NCC, group = group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  topTags(qlf)
  dim(qlf)

  write.csv(qlf$table, file = "SSvsHC_FC.csv")
  
  eset <- log(qlf$fitted.values+1, 2)
  eset[1:3, 1:3]
  write.csv(eset, file = "log2Eset.csv")
}

eset <- read.csv(file = "d:/HP_rdata/HP_treated_PBMC.csv")
qlf <- read.csv(file = "d:/HP_rdata/SSvsHC_FC.csv", row.names = 1)
str(eset)
min(eset$GD_HC.pval)
eset_fc <- eset[eset$GD_HC.fc > 0 & eset$GD_HC.pval < 0.25,]
dim(eset_fc)
#IF_genes <- read.csv(file = "d:/GO_term/cytokine_chemokine2.csv")
IF_genes <- read.csv(file = "d:/GO_term/Hsa_Inflammation_genes.csv")
eset_fc2 <- eset_fc[eset_fc$Gene.symbol %in% c(IF_genes$GOBP_Hsa_Inflammatory.response),]
rownames(qlf)
qlf2 <- qlf[rownames(qlf) %in% eset_fc2$Gene.symbol, ]
qlf3 <- qlf2[order(-qlf2$logFC),]
qlf4 <- qlf3[qlf3$PValue < 0.05, ]

dim(eset_fc2)
str(eset_fc2)
rownames(eset_fc2) <- eset_fc2$Gene.symbol
hm_table <- eset_fc2[,c(13,14,16,17)]
hm_table2 <- hm_table[rownames(qlf4),]

# FUT check
{
  FTF <- c(paste0(rep("FUT", 10), seq(1, 11, 1)), "POFUT1", "POFUT2")
  FTF2 <- c(paste0(rep("FUT", 9), seq(3, 11, 1)), "HP")
  eset_fut <- eset[rownames(eset) %in% FTF2,]
  rownames(eset_fut) <- eset_fut$Gene.symbol
  dim(eset_fut)
  eset_fut[-4,]
  hm_table_fut <- eset_fut[c(1,4,5,6),]
}

# Volcano plot 그리기
{
  qlf1 <- qlf$table
  colnames(qlf1) <- c("Gene", "logFC", "PValue")
  
  count <- dim(qlf1)[1]
  
  for(i in 1:count) {
    if((qlf1$logFC[i] > 1 & qlf1$PValue[i] < 0.05)) {
      qlf1$label[i] <- rownames(qlf1)[i]
    } else {
      qlf1$label[i] <- NA
    }
  }
  
  qlf1$IF_label <- NA
  qlf1$IF_label[qlf1$label %in% c(IF_genes$Cytokine, 
                         IF_genes$Chemokine, 
                         IF_genes$Chemokine.Receptor, "TNF")] <- qlf1[qlf1$label %in% c(IF_genes$Cytokine, 
                                          IF_genes$Chemokine, 
                                          IF_genes$Chemokine.Receptor, "TNF"), 5]
  
  factor(qlf1$IF_label)

  #---SSvsHC---#
  {
    fname = "SSvsHC.pdf"
    
    qlf1$DE <- "ns"
    qlf1$DE[qlf1$logFC > 1 & qlf1$PValue < 0.05] <- "UP"
    qlf1$DE[qlf1$logFC < -1 & qlf1$PValue < 0.05] <- "DOWN"
    qlf1$logpval <- -log10(qlf1$PValue)
    
    Upgenes <- qlf1[qlf1$DE == "UP", ]
    Downgenes <- qlf1[qlf1$DE == "DOWN", ]
    
    library(extrafont)
    
    p <- ggplot(data=qlf1, aes(x=logFC, y=-log10(PValue), fill=DE, label = IF_label)) +
      geom_point(shape = 21, fill = 'white', stroke = 1, aes(color = DE), size = 2) +
      #geom_text(col = "black", hjust=-.1, size = 3, check_overlap=TRUE) +
      geom_vline(xintercept=c(-1), col="black", linetype = 2, size = 0.5) + 
      geom_vline(xintercept=c(1), col="black", linetype = 2, size = 0.5) +
      geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 0.5) +
      scale_color_manual(values=c("#1E90FF", "gray", "#CD5C5C")) +
      scale_x_continuous(limits=c(-12.5,12.5), breaks=c(seq(-11, -1, 5),seq(1, 11, 5))) +
      scale_y_continuous(limits=c(0,10), breaks=c(seq(0, 10, 2.5))) +
      theme_classic() +
      theme(plot.title = element_text(family = "Arial Narrow", hjust = 0.5, size = 30, color = "black")) + # theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))   # 글씨체, 글씨 모양, 가운데 정렬, 크기, 색상을 설정합니다.
      theme(legend.position = "none") +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
      theme(text = element_text(size = 16), #6)
            axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
            axis.text.y = element_text(size = 16, color = "black"), #8)
            axis.line = element_line(size = 0.5, linetype = "solid", color = "black")) + #9)
      xlab(bquote("log"["2"]~"F.C")) +
      ylab(bquote("-log"["10"]~"p value")) +
      geom_label_repel(size = 3.5, colour="black", max.overlaps = getOption("ggrepel.max.overlaps", default = 30), 
                       box.padding = 1) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    #geom_label_repel(size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
    
    ggsave(filename = fname, dpi = 72, width = 375, height = 300, units = 'px') #12)
  }
}

### 히트맵 그리기 ###
{
  cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                  brewer.pal(n = 9, name = "Greys")[3], 
                                  brewer.pal(n = 9, name = "Blues")[7]))(19))
  dev.off()
  
  a <- pheatmap(hm_table2, 
                scale = "row", 
                legend = T, 
                color = cols,
                border_color = "black", 
                #gaps_col = c(12, 31),
                clustering_distance_rows = "euclidean",
                fontsize = 7,
                #gaps_row = 1:dim(mf7)[1],
                cellwidth = 15, 
                cellheight = 10,
                legend_breaks = c(-1, 0, 1),
                breaks = seq(-1, 1, 2/19),
                show_colnames = T,
                show_rownames = T,
                cluster_rows = F,
                cluster_cols = F, 
                angle_col = 0, 
                family = "Arial Narrow")
}

# GO analysis 그리기
{
  organism <- "org.Hs.eg.db"
  library(organism, character.only = TRUE)
  library(clusterProfiler)
  library(enrichplot)
  
  ### G vs HC ###
  {
    
    eset_fc3 <- eset_fc[order(-eset_fc$GD_HC.fc),]
    eset_fc4 <- eset_fc3[!duplicated(eset_fc3$Gene.symbol),]
    rownames(eset_fc4) <- eset_fc4$Gene.symbol
    geneList <- eset_fc4
    dim(geneList)
    
    original_gene_list <- geneList$GD_HC.fc
    names(original_gene_list) <- rownames(geneList)
    length(original_gene_list)
    gene_list <- na.omit(original_gene_list)
    length(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
    dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
    
    df2 <- gene_list
    df3 <-ids[ids$SYMBOL %in% names(df2), 2]
    names(df2) <- df3
    length(df2)
    df2 <- sort(df2, decreasing = TRUE)
    
    gse <- gseGO(geneList     = df2,
                 ont = "BP", 
                 organism,
                 keyType       = "ENTREZID", 
                 nPerm        = 10000,
                 minGSSize    = 10,
                 maxGSSize    = 60,
                 pvalueCutoff = 0.99,
                 pAdjustMethod = "none",
    )
    
    final <- as.data.frame(gse@result)
    str(final)
    final$`'-log2(adj.pvalue)` <- -log(final$p.adjust, 2)
    
    require(DOSE)
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
    dotplot(gse, showCategory=20, split=".sign", label_format = 50, decreasing = F
    ) + facet_grid(.~.sign)
    
    path <- c("response to chemokine", "cellular response to chemokine", "chemokine-mediated signaling pathway", 
              "response to interleukin-1", "cellular response to interleukin-1", 'regulation of leukocyte migration', 
              "neutrophil chemotaxis", "neutrophil migration")
    
    library(enrichplot)
    dotplot(gse, x = "GeneRatio", color = 'p.adjust', showCategory=path, label_format = 50) +
    scale_x_continuous(limits = c(0, 0.7))
  }
}
