#install.packages('extrafont')
library(extrafont)
font_import(pattern = "Arial Narrow")

install.packages("ggsignif")                      # Install ggsignif package
library("ggsignif")   

temp <- getwd()
setwd(temp)
setwd("C:/Users/User/Desktop/01_Sepsis_Hp_Project/GSE65682/231228")
qlf1 <- read.csv(file = "GSE65682_FC_table.csv", row.names = 1)
qlf2 <- qlf1[order(-qlf1$logFC),]
qlf3 <- qlf2[!(duplicated(qlf2$gene)), ]
qlf4 <- qlf3[(qlf3$logFC > 0 & qlf3$PValue < 0.05),]
dim(qlf4)
eset <- read.csv(file = "GSE65682_Eset.csv", row.names = 1)
eset[1:3, 1:3]
str(eset)
dim(eset)
eset2 <- eset[rownames(eset) %in% rownames(qlf4),]
eset2[1:3, 1:3]
dim(eset2)
rownames(eset2)
Gene <- qlf4[rownames(eset2),"gene"]
eset3 <- cbind(Gene, eset2)
eset3[1:10, 1:3]
eset4 <- eset3[,-1]
rownames(eset4) <- eset3$Gene
eset4[1:10, 1:3]
eset5 <- as.data.frame(t(eset4))
str(eset4)
eset_HC <- eset5[1:42,]
eset_SS <- eset5[43:dim(eset5)[1],]

# Box plot 그리기

  {
    target <- "HP"
    
    MinY <- floor(min(eset5[,target]))
    MaxY <- ceiling(max(eset5[,target])) + 2
    
    A <- max(eset_HC[,target])
    B <- max(eset_SS[,target])

    if(A > B) {
      YP1 <- A + 0.3
    } else {
      YP1 <- B + 0.3
    } 
    FYP <- c(YP1)
  }
  
temp_expr <- c(eset_HC[,target], eset_SS[,target])
length(temp_expr)
temp_group <- c(rep("HC",length(eset_HC[,target])), rep("SS",length(eset_SS[,target])))
length(temp_group)
temp_table <- as.data.frame(cbind(Group = temp_group, Expr = temp_expr))
temp_table$Expr <- as.numeric(temp_table$Expr)
temp_table2 <- temp_table[!(temp_table$Expr == 0),]

temp_pval <- qlf4[qlf4$gene == target, "PValue"]

  if(temp_pval < 0.001) {
    p1 <- "***"
  } else if(temp_pval < 0.01) {
    p1 <- "**"
  } else if(temp_pval < 0.05) {
    p1 <- "*" 
  } else {
    p1 <- "ns"
  }
  
  anot <- p1
  
  p <- ggplot(temp_table, aes(y=`Expr`, x=`Group`, fill=`Group`)) +  
    geom_boxplot(position=position_dodge(1), color = "black", linewidth = 0.5)
  
  # 그래프에 1) 점 표시 / 2) X 축 라벨 지정 / 3) Y 축 범위 및 간격 지정 / 4) 그래프 타이틀 및  X축, Y축 이름 지정 / 
  # 5) 그래프 색지정 / 6) 레전드 제거 및 글씨체 및 글씨크기 설정 / 7), 8) xy축 글씨 크기 및 색깔 위치 지정 / 
  # 9) 축 라인 굵기, 스타일, 색 설정 / 10) 축의 틱 굵기, 길이, 색 설정 / 11) 유의수준 표시 / 12) 그래프 저장
  
  library(ggsignif)
  
  p + 
    geom_dotplot(binaxis='y', stackdir='center',dotsize = 0, position=position_dodge(1), color = "black", fill = "black") + #1)
    scale_x_discrete(limits=c("HC", "SS")) + #2)
    #scale_y_continuous(limits = c(FminY, FmaxY), breaks = seq(FminY, FmaxY, tickgap)) + #3)
    labs(title = target, x="", y = "Expression level\n", face = "italic") + #4)
    scale_fill_grey() + theme_classic() + #5)
    theme(legend.position="none", text = element_text(size = 16, family="Arial Narrow"), #6)
          axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
          axis.text.y = element_text(size = 16, color = "black"), #8)
          axis.line = element_line(size = 0.5, linetype = "solid", color = "black"), #9)
          axis.ticks = element_line(size = 0.5, linewidth = 10, linetype = "solid", color = "black"), 
          plot.title = element_text(face = "italic")) + #10)
    geom_signif(comparisons = list(c("HC", "SS")), #11)
                map_signif_level = TRUE, 
                y_position = FYP,
                size = 0.5, 
                family = "Arial Narrow",
                textsize = 5, 
                annotations = anot)

  ggsave(filename = paste0(target, ".pdf"), dpi = 72, width = 150, height = 300, units = 'px') #12)
 

# Volcano plot 그리기
{

  
  #---SSvsHC---#
  {
    library(ggplot2)
    
    fname = "SSvsHC.pdf"
    
    qlf1$DE <- "ns"
    qlf1$DE[qlf1$logFC > 1 & qlf1$PValue < 0.05] <- "UP"
    qlf1$DE[qlf1$logFC < -1 & qlf1$PValue < 0.05] <- "DOWN"
    
    Upgenes <- qlf1[qlf1$DE == "UP", ]
    Downgenes <- qlf1[qlf1$DE == "DOWN", ]
    
    p <- ggplot(data=qlf1, aes(x=logFC, y=-log10(PValue), col=DE)) +
      geom_point(size = 3) + 
      #geom_text(col = "black", hjust=-.1, size = 3, check_overlap=TRUE) +
      geom_vline(xintercept=c(-1), col="black", linetype = 2, size = 0.5) + 
      geom_vline(xintercept=c(1), col="black", linetype = 2, size = 0.5) +
      geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 0.5) +
      scale_color_manual(values=c("#1E90FF", "gray", "#CD5C5C")) +
      scale_x_continuous(limits=c(-12.5,12.5), breaks=c(seq(-11, -1, 5),seq(1, 11, 5))) +
      scale_y_continuous(limits=c(0,15), breaks=c(seq(0, 15, 2.5))) +
      theme_classic() +
      theme(plot.title = element_text(family = "Arial Narrow", hjust = 0.5, size = 30, color = "black")) + # theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))   # 글씨체, 글씨 모양, 가운데 정렬, 크기, 색상을 설정합니다.
      theme(legend.position = "none") +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
      theme(text = element_text(size = 16), #6)
            axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
            axis.text.y = element_text(size = 16, color = "black"), #8)
            axis.line = element_line(size = 0.5, linetype = "solid", color = "black")) + #9)
      xlab(bquote("log"["2"]~"F.C")) +
      ylab(bquote("-log"["10"]~"p value"))
    #geom_label_repel(size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
    
    ggsave(filename = fname, dpi = 72, width = 375, height = 300, units = 'px') #12)
  }

  
}

# GO analysis 그리기
{
  organism <- "org.Hs.eg.db"
  library(organism, character.only = TRUE)
  library(clusterProfiler)
  library(enrichplot)
  
  ### SS vs HC ###
  {
    geneList <- qlf1[qlf1$PValue<0.05,]
    dim(geneList)
    
    original_gene_list <- geneList$logFC
    names(original_gene_list) <- rownames(geneList)
    length(original_gene_list)
    gene_list <- na.omit(original_gene_list)
    length(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
    dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
    
    df2 = gene_list[dedup_ids$SYMBOL]
    df3 <-ids[ids$SYMBOL %in% names(df2), 2]
    
    names(df2) <- df3
    length(df2)
    df2 <- sort(df2, decreasing = TRUE)
    
    gse <- gseKEGG(geneList     = df2,
                   organism     = "hsa",
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 100,
                   pvalueCutoff = 0.99,
                   pAdjustMethod = "none",
                   keyType       = "kegg")
    
    final <- as.data.frame(gse@result)
    str(final)
    final$`'-log2(adj.pvalue)` <- -log(final$p.adjust, 2)
    
    require(DOSE)
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
    dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
    
    enrich <- strsplit(final$core_enrichment, split="/")
    
    for(i in 1:length(enrich)) {
      if(i == 1) {
        temp <- length(enrich[[i]])
      }
      else {
        temp <- c(temp, length(enrich[[i]]))
      }
    }
    
    length(enrich[[1]])
    
    final$count <- temp
    write.csv(final, file = "SSvsH_seq_KEGG_result.csv")
  }
  #########################################

  
}

# GSE65682 분석
### edgeR로 DEG 구하기 ###
{
  #BiocManager::install('edgeR')
  
  library(edgeR)
  
  NCC <- readRDS(file = "d:/GSE65682/GSE65682_eset.rds")
  HC <- readRDS(file = "d:/GSE65682/HC_eset.rds")  ### 
  G <- readRDS(file = "d:/GSE65682/SP_eset.rds")
  B <- readRDS(file = "d:/GSE65682/SK_eset.rds")
  gene <- readRDS(file = "d:/GSE65682/Gene_annotation.rds")
  
  length(colnames(HC))
  length(colnames(G))
  length(colnames(B))
  
  colnames(NCC)[which(colnames(NCC) %in% colnames(HC))] <- paste0(rep("HC", 42), sprintf("%03d", seq(1,42,1)))
  colnames(NCC)[which(colnames(NCC) %in% colnames(G))] <- paste0(rep("G", 365), sprintf("%03d", seq(1,365,1)))
  colnames(NCC)[which(colnames(NCC) %in% colnames(B))] <- paste0(rep("B", 114), sprintf("%03d", seq(1,114,1)))
  View(NCC)

  group <- factor(c(rep('HC', 42), rep('SS', 479)))
  y <- DGEList(counts = NCC, group = group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  design <- model.matrix(~group)
  colnames(design) <- levels(y$samples$group)
  design
  
  y <- estimateDisp(y)
  
  fit <- glmQLFit(y, design)
  qlf2 <- glmQLFTest(fit,coef=2)
  topTags(qlf2)

  write.csv(qlf2$table, file = "GSE65682_FC_SSvsHC.csv")
  write.csv(y$counts, file = "GSE65682_Eset.csv")
  
  fc_table <- qlf2$table
  
  dim(fc_table)
  
  fc_table$DE <- "ns"
  fc_table$DE[fc_table$logFC > 0 & fc_table$PValue < 0.05] <- "UP"
  fc_table$DE[fc_table$logFC < 0 & fc_table$PValue < 0.05] <- "DOWN"
  
  fc_table$gene <- gene[gene$V2 %in% rownames(fc_table), 1]
  
  UP2 <- fc_table[fc_table$DE == "UP",]
  DOWN2 <- fc_table[fc_table$DE == "DOWN",]
  
  UP2_gene <- gene[gene$V2 %in% rownames(UP2), 1]
  DOWN2_gene <- gene[gene$V2 %in% rownames(DOWN2), 1]
  
  length(UP2_gene)
  length(DOWN2_gene)
  DOWN2_gene[120:135] <- NA
  deg_table <- cbind(UP2_gene, DOWN2_gene)
  write.csv(deg_table, file = "GSE65682_UPDOWN_Genes.csv")
  write.csv(fc_table, file = "GSE65682_FC_table.csv")
}


# GSE65682 분석
### Pathway 분석 ###
{
  organism <- "org.Hs.eg.db"
  library(organism, character.only = TRUE)
  library(clusterProfiler)
  library(enrichplot)
  
  ### SS vs HC ###
  {
    qlf2 <- fc_table
    geneList <- qlf2[qlf2$PValue<0.05,]
    dim(geneList)
    
    original_gene_list <- geneList$logFC
    names(original_gene_list) <- geneList$gene
    length(original_gene_list)
    gene_list <- na.omit(original_gene_list)
    length(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
    dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
    dim(dedup_ids)
    which(duplicated(dedup_ids))
    ids$SYMBOL[which(duplicated(ids$SYMBOL))]
    which((dedup_ids$SYMBOL %in% ids$SYMBOL))
    df2 = gene_list[dedup_ids$SYMBOL]
    
    ids <- ids[!duplicated(ids$SYMBOL), ]
    
    df3 <- ids[ids$SYMBOL %in% names(df2), 2]
    length(df2)
    length(df3)
    
    names(df2) <- df3
    df2 <- sort(df2, decreasing = TRUE)
    
    gse <- gseKEGG(geneList     = df2,
                   organism     = "hsa",
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 100,
                   pvalueCutoff = 0.99,
                   pAdjustMethod = "none",
                   keyType       = "kegg")
    
    final <- as.data.frame(gse@result)
    str(final)
    final$`'-log2(adj.pvalue)` <- -log(final$p.adjust, 2)
    
    require(DOSE)
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
    dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
    
    enrich <- strsplit(final$core_enrichment, split="/")
    
    for(i in 1:length(enrich)) {
      if(i == 1) {
        temp <- length(enrich[[i]])
      }
      else {
        temp <- c(temp, length(enrich[[i]]))
      }
    }
    
    length(enrich[[1]])
    
    final$count <- temp
    write.csv(final, file = "GvsH_seq_KEGG_result.csv")
  }
  #########################################
  
  ### B vs HC ###
  {
    geneList <- qlf2[qlf2$PValue<0.05,]
    dim(geneList)
    
    original_gene_list <- geneList$logFC
    names(original_gene_list) <- rownames(geneList)
    length(original_gene_list)
    gene_list <- na.omit(original_gene_list)
    length(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
    dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
    
    df2 = gene_list[dedup_ids$SYMBOL]
    df3 <-ids[ids$SYMBOL %in% names(df2), 2]
    
    names(df2) <- df3
    length(df2)
    df2 <- sort(df2, decreasing = TRUE)
    
    gse <- gseKEGG(geneList     = df2,
                   organism     = "hsa",
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 100,
                   pvalueCutoff = 0.99,
                   pAdjustMethod = "none",
                   keyType       = "kegg")
    
    final <- as.data.frame(gse@result)
    str(final)
    final$`'-log2(adj.pvalue)` <- -log(final$p.adjust, 2)
    
    require(DOSE)
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
    dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
    
    enrich <- strsplit(final$core_enrichment, split="/")
    
    for(i in 1:length(enrich)) {
      if(i == 1) {
        temp <- length(enrich[[i]])
      }
      else {
        temp <- c(temp, length(enrich[[i]]))
      }
    }
    
    length(enrich[[1]])
    
    final$count <- temp
    write.csv(final, file = "BvsHC_seq_KEGG_result.csv")
  }
  #################################################
  
  
  ### G vs B ###
  {
    geneList <- qlf3[qlf3$PValue<0.05,]
    dim(geneList)
    
    original_gene_list <- geneList$logFC
    names(original_gene_list) <- rownames(geneList)
    length(original_gene_list)
    gene_list <- na.omit(original_gene_list)
    length(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
    dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
    
    df2 = gene_list[dedup_ids$SYMBOL]
    df3 <-ids[ids$SYMBOL %in% names(df2), 2]
    
    names(df2) <- df3
    length(df2)
    df2 <- sort(df2, decreasing = TRUE)
    
    gse <- gseKEGG(geneList     = df2,
                   organism     = "hsa",
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 100,
                   pvalueCutoff = 0.99,
                   pAdjustMethod = "none",
                   keyType       = "kegg")
    
    final <- as.data.frame(gse@result)
    str(final)
    final$`'-log2(adj.pvalue)` <- -log(final$p.adjust, 2)
    
    require(DOSE)
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
    dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
    
    enrich <- strsplit(final$core_enrichment, split="/")
    
    for(i in 1:length(enrich)) {
      if(i == 1) {
        temp <- length(enrich[[i]])
      }
      else {
        temp <- c(temp, length(enrich[[i]]))
      }
    }
    
    length(enrich[[1]])
    
    final$count <- temp
    write.csv(final, file = "GvsB_seq_KEGG_result.csv")
  }
  #################################################
  
  
}



# Venn diagram 그리기 (다른 2개 코호트)
{
  qlf2 <- 

  
  #---GvsHC---#
  {
    qlf1$DE <- "ns"
    qlf1$DE[qlf1$logFC > 1 & qlf1$PValue < 0.05] <- "UP"
    qlf1$DE[qlf1$logFC < -1 & qlf1$PValue < 0.05] <- "DOWN"
    
    Upgenes1 <- qlf1[qlf1$DE == "UP", ]
    Downgenes1 <- qlf1[qlf1$DE == "DOWN", ]
  }
  
  #---BvsHC---#
  {
    qlf2$DE <- "ns"
    qlf2$DE[qlf2$logFC > 1 & qlf2$PValue < 0.05] <- "UP"
    qlf2$DE[qlf2$logFC < -1 & qlf2$PValue < 0.05] <- "DOWN"
    
    Upgenes2 <- qlf2[qlf2$DE == "UP", ]
    Downgenes2 <- qlf2[qlf2$DE == "DOWN", ]
  }
  
  #---BvsG---#
  {
    qlf3$DE <- "ns"
    qlf3$DE[qlf3$logFC > 1 & qlf3$PValue < 0.05] <- "UP"
    qlf3$DE[qlf3$logFC < -1 & qlf3$PValue < 0.05] <- "DOWN"
    
    Upgenes3 <- qlf3[qlf3$DE == "UP", ]
    Downgenes3 <- qlf3[qlf3$DE == "DOWN", ]
  }
  
  #---Venn diagram---#
  {
    library(VennDiagram)
    library(RColorBrewer)
    
    coln <- length(x)
    
    x = list(
      `GvsHC` = rownames(Upgenes1),
      `BvsHC` = rownames(Upgenes2), 
      `GvsB` = rownames(Downgenes3) 
    )
    
    Up_overlap <- calculate.overlap (x)
    
    write.csv(Up_overlap$a5, file = "Up_overlap.csv")
    saveRDS(Up_overlap, file = "UP_overlap.rds")
    
    a <- venn.diagram(x,
                      filename = 'VD_UP.tiff',   # 저장할 파일명을 설정합니다.
                      col = "white",   # 벤 다이어그램 테두리 색상을 설정합니다.
                      lty = 1,   # 벤 다이어그램 테두리 선 모양을 설정합니다.
                      fill = c(brewer.pal(n = 9, name = 'Set1')[1:coln]),   # 각 벤 다이어그램 내부 색상을 설정합니다.
                      alpha = 0.20,   # 벤 다이어그램의 알파 투명도를 설정합니다.
                      main.fontfamily = "Arial Narrow",
                      cat.fontfamily = "Arial Narrow",
                      cat.col = "black",   # 각 벤 다이어그램의 명칭에 대한 색상을 설정합니다.
                      cat.cex = 1,   # 각 벤 다이어그램의 명칭에 대한 글자 크기를 설정합니다.
                      cat.fontface = "bold",   # 각 벤 다이어그램의 명칭에 대한 글자를 볼드체로 설정합니다.
                      margin = 0.1,   #  벤 다이어그램 주위의 공간을 설정합니다.
                      scaled = F)
    
  }
  
  x = list(
    `GvsHC` = rownames(Downgenes1),
    `BvsHC` = rownames(Downgenes2), 
    `GvsB` = rownames(Upgenes3) 
  )
  
  coln <- length(x)
  
  Down_overlap <- calculate.overlap (x)
  
  write.csv(Down_overlap$a5, file = "Down_overlap.csv")
  saveRDS(Down_overlap, file = "DOWN_overlap.rds")
  
  a <- venn.diagram(x,
                    filename = 'VD_DOWN.tiff',   # 저장할 파일명을 설정합니다.
                    col = "white",   # 벤 다이어그램 테두리 색상을 설정합니다.
                    lty = 1,   # 벤 다이어그램 테두리 선 모양을 설정합니다.
                    fill = c(brewer.pal(n = 9, name = 'Set1')[1:coln]),   # 각 벤 다이어그램 내부 색상을 설정합니다.
                    alpha = 0.20,   # 벤 다이어그램의 알파 투명도를 설정합니다.
                    main.fontfamily = "Arial Narrow",
                    cat.fontfamily = "Arial Narrow",
                    cat.col = "black",   # 각 벤 다이어그램의 명칭에 대한 색상을 설정합니다.
                    cat.cex = 1,   # 각 벤 다이어그램의 명칭에 대한 글자 크기를 설정합니다.
                    cat.fontface = "bold",   # 각 벤 다이어그램의 명칭에 대한 글자를 볼드체로 설정합니다.
                    margin = 0.1,   #  벤 다이어그램 주위의 공간을 설정합니다.
                    scaled = F)
  
}
