### GSE232753 ###

library(edgeR)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE232753", "file=GSE232753_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID
annot[1:3, 1:3]

library(edgeR)
### edgeR로 DEG 구하기 ###
{
  group <- factor(c(rep('HC', 8), rep('SS', 20)))
  y <- DGEList(counts = tbl, group = group)
  design <- model.matrix(~group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  topTags(qlf)
  dim(qlf)
  
  gene_symbol <- annot[rownames(annot) %in% rownames(qlf),c(1,2)]
  head(rownames(qlf))
  head(rownames(gene_symbol))
  
  dim(qlf$table)
  fc_table <- qlf$table
  fc_table2 <- cbind(Gene = gene_symbol, fc_table)
  write.csv(fc_table2, file = "GSE232753_SSvsHC.csv")
  
  log2eset <- log(qlf$fitted.values+1, 2)
  boxplot(log(eset+1, 2))
  colnames(log2eset) <- c(paste0(rep("HC", 8), sprintf("%02d", seq(1,8,1))), 
                          paste0(rep("SS", 20), sprintf("%02d", seq(1,20,1))))
  eset <- qlf$fitted.values
  genetable <- annot[rownames(annot) %in% rownames(eset),c(1,2)]
  eset2 <- cbind(genetable, eset)
  rownames(eset2) <- eset2$Symbol
  dim(eset2)
  eset2[1:3, 1:3]
  eset3 <- eset2[,3:dim(eset2)[2]]
  colnames(eset3) <- c(paste0(rep("HC", 8), sprintf("%02d", seq(1,8,1))), 
                       paste0(rep("SS", 20), sprintf("%02d", seq(1,20,1))))
  
  
  
  gene_symbol2 <- annot[rownames(annot) %in% rownames(log2eset),c(1,2)]
  
  dim(log2eset)
  dim(gene_symbol2)
  
  log2eset2 <- cbind(gene_symbol2, log2eset)
  write.csv(log2eset2, file = "GSE232753_log2eset.csv")
}  

### PCA plot ###
{
  colnames(y$counts) <- c(paste0(rep("HC", 8), sprintf("%02d", seq(1,8,1))), 
                           paste0(rep("SS", 20), sprintf("%02d", seq(1,20,1))))
  shape <- c(rep(1, 8), rep(1, 20))
  plotMDS(y, col=c(rep("black",8), rep("red",20)))
  plotMDS(y, col=c(rep("red",24), rep("black",12)), pch=shape)
  legend("topright", col=c("black","red"), pch=c(1,1), legend=c("HC","Sepsis"))
  
  getwd()
}

### insterested genes ###
target_genes <- c("IL1B", "IL6", "TNF", "HP", "CLEC4E")
target_genes2 <- c(paste0(rep("CXCL", 30), seq(1,30,1)), paste0(rep("CCL", 30), seq(1,30,1)), "S100A8", "S100A9")
FTF <- c(paste0(rep("FUT", 10), seq(1, 11, 1)), "POFUT1", "POFUT2")
NAT <- c("A4GNT", "ALG13", "ALG14", "DPAGT1", "EOGT", "GNPTAB", "GNPTG", 
         "MGAT2", "MAGAT3", "MGAT4A", "MGAT4B", "MGAT5", "OGT")
final_target_genes <- c(target_genes, target_genes2, FTF, NAT)
length(final_target_genes)

target_eset <- log2eset2[log2eset2$Symbol %in% final_target_genes,]
rownames(target_eset) <- target_eset$Symbol
target_eset <- target_eset[,-1]
target_eset <- target_eset[,-1]
dim(target_eset)

fc_table3 <- fc_table2[fc_table2$Gene.Symbol %in% final_target_genes,]
rownames(fc_table3) <- fc_table3$Gene.Symbol

# WT_SS Interested gene graph 그리기
for(i in 1:dim(target_eset)[1]) {
  
  fname1 <- paste0(sprintf("%02d", i), "_WTvsSS_", rownames(target_eset)[i], ".pdf")
  title <- rownames(target_eset)
  
  temp_expr <- target_eset[i,]
  length(temp_expr)
  temp_group <- c(rep("WT",8), rep("SS",20))
  length(temp_group)
  temp_table <- as.data.frame(cbind(Group = temp_group, Expr = as.vector(temp_expr)))
  temp_table$Expr <- as.numeric(temp_table$Expr)
  temp_table$Group <- temp_group
  temp_table$Group = factor(x = temp_table$Group, levels = c("WT", "SS"))
  
  # p value 위치설정
  
  {
    WT_pos3 <- target_eset[,1:8]
    KO_pos3 <- target_eset[,9:28]
    A <- max(WT_pos3[i,])
    B <- max(KO_pos3[i,])
    
    MinY <- floor(min(temp_table[,2]))
    MaxY <- ceiling(max(temp_table[,2])) + 2
    
    gap <- (MaxY - MinY) / 5
    
    if(A > B) {
      YP1 <- A + 0.3
    } else {
      YP1 <- B + 0.3
    } 
    FYP <- c(YP1)
  }
  
  temp_pval <- fc_table3$PValue[i]
  
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
  
  
  p + 
    geom_dotplot(binaxis='y', stackdir='center',dotsize = 0, position=position_dodge(1), color = "black", fill = "black") + #1)
    scale_x_discrete(limits=c("WT", "SS")) + #2)
    #scale_y_continuous(limits = c(FminY, FmaxY), breaks = seq(FminY, FmaxY, tickgap)) + #3)
    labs(title = title[i], x="", y = "Expression level\n", face = "italic") + #4)
    scale_fill_grey() + theme_classic() + #5)
    theme(legend.position="none", text = element_text(size = 16, family="Arial Narrow"), #6)
          axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
          axis.text.y = element_text(size = 16, color = "black"), #8)
          axis.line = element_line(size = 0.5, linetype = "solid", color = "black"), #9)
          axis.ticks = element_line(size = 0.5, linewidth = 10, linetype = "solid", color = "black"), 
          plot.title = element_text(face = "italic")) + #10)
    geom_signif(comparisons = list(c("WT", "SS")), #11)
                map_signif_level = TRUE, 
                y_position = FYP,
                size = 0.5, 
                family = "Arial Narrow",
                textsize = 5, 
                annotations = anot)
  
  
  ggsave(filename = fname1, dpi = 72, width = 150, height = 300, units = 'px') #12)
  
}

#########################


log2eset3 <- as.data.frame(t(eset3))
log2eset3[1:3,1:3]
HC_leset <- log2eset3[1:8,]
dim(HC_leset)
SS_leset <- log2eset3[9:28,]
dim(SS_leset)

basegene <- "CLEC4E"
final_target_genes <- c("LYZ", "IL10", "ACTINB")

#######################

{
  for(i in 1:dim(SS_leset)[2]) {
    if(i == 1) {
      temp <- cor.test(SS_leset[,basegene], SS_leset[,i])
      ctable <- as.data.frame(cbind(corr = temp$estimate, p.value = temp$p.value))
      rownames(ctable) <- colnames(SS_leset)[i]
      
    } else {
      temp <- cor.test(SS_leset[,basegene], SS_leset[,i])
      ttable <- as.data.frame(cbind(corr = temp$estimate, p.value = temp$p.value))
      ctable <- as.data.frame(rbind(ctable, ttable))
      rownames(ctable)[i] <- colnames(SS_leset)[i]
    }
  }
  
  write.csv(ctable, file = "SS_CLEC4E_Corr.csv")
  saveRDS(ctable, file = "SS_CLEC4E_Corr.rds")
}  

# SS Correlation
final_target_genes <- c("LYZ", "IL10", "ACTINB")
control <- "CLEC4E"
count <- 0
for(i in final_target_genes) {
  
  if(!(control %in% log2eset2$Symbol)) {
    next
  } 
  
  if(!(i %in% log2eset2$Symbol)) {
    next
  } 
  
  count <- count + 1
  con_expr <- as.numeric(log2eset2[log2eset2$Symbol == control,3:30])
  tar_expr <- as.numeric(log2eset2[log2eset2$Symbol == i,3:30])

  fname2 <- paste0(sprintf("%02d", count), "_KO_Corr_", i, ".pdf")
  
  SScon <- con_expr[9:28]
  SSset <- tar_expr[9:28]
  
  temp <- cor.test(SScon, SSset)
  
  EP <- cbind(Rho = temp$estimate, `p value` = temp$p.value)
  rownames(EP) <- i
  
  if(i == head(final_target_genes,1)){
    data <- EP 
  } else {
    data <- rbind(data, EP)
  }
  
  if(round(data[i,2], 4) == 0) {
    pval <- "0.000"
  } else {
    pval <- round(data[i,2], 4)
  }
  
  rlabel <- paste0("Rho = ", signif(round(data[i,1], 4)), '\n', "p = ", pval, " ")
  
  temp_corr <- as.data.frame(cbind(SScon, SSset))
  
  ggplot(temp_corr, aes(x = SScon, y = SSset)) +
    geom_point(col = "black", size=2) +
    theme_bw() +
    xlab(control) +
    ylab(i) +
    geom_smooth(color = "#0066CC", method = "lm") +
    annotate("text", colour="red", size=5, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2, family="Arial Narrow") +
    theme(legend.position="none", text = element_text(size = 16, family="Arial Narrow"), #6)
          axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
          axis.text.y = element_text(size = 16, color = "black"), #8)
          axis.line = element_line(size = 0.5, linetype = "solid", color = "black"), #9)
          axis.ticks = element_line(size = 0.5, linewidth = 10, linetype = "solid", color = "black"), 
          plot.title = element_text(face = "italic"))
  
  ggsave(filename = fname2, dpi = 72, width = 300, height = 300, units = 'px') #12)
  
  if(i == tail(final_target_genes, 1)) {
    
    write.csv(data, file = "GSE232753_Correlation.csv")
  }
}
