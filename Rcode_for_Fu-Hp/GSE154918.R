G154 <- readRDS(file = "d:/GSE154918_eset.rds")
head(G154)
dim(G154)
colnames(G154)

G154_HC <- readRDS(file = "d:/gse154918/HC_eset.rds")
G154_SP <- readRDS(file = "d:/gse154918/SP_eset.rds")
G154_SP2 <- readRDS(file = "d:/gse154918/SK_eset.rds")

dim(G154_HC)
dim(G154_SP)
dim(G154_SP2)

G154 <- cbind(G154_HC, G154_SP, G154_SP2)
G154[1:5, 1:5]
colnames(G154) <- c(paste0(rep("HC", 40), sprintf("%02d", seq(1,40,1))), paste0(rep("SS", 39), sprintf("%02d", seq(1,39,1))))

### edgeR로 DEG 구하기 ###
{
  library(edgeR)
  
  newG154 <- cbind(G154[,41:dim(G154)[2]], G154[,1:40])
  group <- factor(c(rep('SS', 39), rep('HC', 40)))
  y <- DGEList(counts = newG154, group = group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  topTags(qlf)
  dim(qlf)
  
  shape <- c(rep(1, 39), rep(1, 40))
  plotMDS(y, col=c(rep("red",39), rep("black",40)))
  plotMDS(y, col=c(rep("red",39), rep("black",40)),pch=shape)
  
  legend("topright", col=c("black","red"), pch=c(1,1), legend=c("HC","Sepsis"))
  
  write.csv(qlf$table, file = "SSvsHC.csv")
  
  eset <- y$counts
  eset[1:3, 1:3]
  write.csv(eset, file = "GSE154918_Eset.csv")
}

#BiocManager::install('gprofiler2')
library(gprofiler2)
ensembl.genes <- rownames(qlf)
gene.table <- gconvert(ensembl.genes,organism="hsapiens",target="ENTREZGENE",filter_na = F)

library(stringr)
gene.table2 <- gene.table[str_detect(gene.table$target_number,'.1'),]
gene.table3 <- gene.table2[!duplicated(gene.table2$input),]

qlf$table$Gene <- gene.table3$name
qlf$table$Gene_ID <- gene.table3$input

write.csv(qlf$table, file = "GSE154918_FC_SSvsHC.csv")

#### Graph ####

target2 <- c("HP", "CLEC4E")

G154_SS <- as.data.frame(t(newG154[,1:39]))
G154_HC <- as.data.frame(t(newG154[,40:dim(newG154)[2]]))
tG154 <- as.data.frame(t(newG154))

rqlf <- qlf$table
tG154[1:3, 1:3]
G154_SS[1:3, 1:3]
G154_HC[,target]

which(colnames(G154_HC) %in% target)

target <- qlf$table$Gene_ID[qlf$table$Gene %in% target2]
target
rqlf[target,4]

#install.packages('extrafont')
library(extrafont)
library(ggplot2)
library(ggsignif)

font_import()
fonts()
font_import(path = "C:/Windows/Fonts", pattern = "Arial Narrow.TTF")
windowsFonts()
extrafont::font_import(pattern="Arial", prompt=FALSE) 

target2 <- c("HP", "CLEC4E")

SPGgenes <- c("S100A8", "FCER1G", "C19orf59", "NFIL3", "H2AFJ", "S100A9", "CLEC4D", "SERPINA1", 
              "S100A12", "FAM20A", "GGH", "GBA", "FLOT2", "HIST1H2BD", "HPR84", "DRAM1", "TKTL1", 
              "AGTRAP", 
              "IRAK3", "BASP1")
MMPgenes <- c("MMP9", "MMP8", "MMP25")
IFgenes <- c("ALOX5", "IL1B", "IL6", "TNF", "IL18", "CCL2", "CCL3", "CCL4", 
             "CXCL8", "CXCL10","CLEC4E", "NLRP3", "CASP1")

FTF <- c(paste0(rep("FUT", 10), seq(1, 11, 1)), "POFUT1", "POFUT2")
NAT <- c("A4GNT", "ALG13", "ALG14", "DPAGT1", "EOGT", "GNPTAB", "GNPTG", 
         "MGAT2", "MAGAT3", "MGAT4A", "MGAT4B", "MGAT5", "OGT")


for(i in c(SPGgenes, MMPgenes, IFgenes, FTF, NAT)) {
  
  if(!(i %in% qlf$table$Gene)) {
    next
  }
  
  target <- qlf$table$Gene_ID[qlf$table$Gene %in% i]

  ### 파일명 설정 ###
  fname1 <- paste0(sprintf("%05d", which(colnames(G154_HC) %in% target)), "_", i, ".pdf")
  
  ### 유의수준 설정 ###
  {
    if(rqlf[target, 4] < 0.001) {
      p1 <- "***"
    } else if(rqlf[target, 4] < 0.01) {
      p1 <- "**"
    } else if(rqlf[target, 4] < 0.05) {
      p1 <- "*" 
    } else {
      p1 <- "ns"
    }
    
    anot <- c(p1)
    
  }
  
  ### p value 위치 설정 ###
  {
    A <- max(G154_SS[,target])
    B <- max(G154_HC[,target])
    
    MinY <- floor(min(tG154[,target]))
    MaxY <- ceiling(max(tG154[,target])) + 2
    
    gap <- (MaxY - MinY) / 5
    
    if(A > B) {
      YP1 <- A + 0.3
      YP2 <- A + gap
    } else {
      YP1 <- B + 0.3
      YP2 <- B + gap
    } 
    
    FYP <- c(YP1, YP2)
  }
  #########################
  
  ### 그래프 y축 설정 ###
  {
    if((MaxY - MinY) / 5 > 1.5 & (MaxY - MinY) / 5 < 2) {
      FmaxY <- MaxY + (2 - (MaxY %% 2))
      FminY <- MinY - (MinY %% 2)
      tickgap <- 2
    } else if((MaxY - MinY) / 5 >= 2) {
      FmaxY <- MaxY + (5 - (MaxY %% 5))
      FminY <- MinY - (MinY %% 5)
      tickgap <- 5
    } else {
      FmaxY <- MaxY
      FminY <- MinY
      tickgap <- 1
    }
  }
  #######################
  
  ### ggplot 형식으로 table 만들기 ###
  {
    temp1 <- c(G154_HC[,target], G154_SS[,target])
    temp2 <- c(rep("HC", 40), rep("SS", 39))
    temp3 <- as.data.frame(cbind(`Expression level` = temp1, `Group` = temp2))
    temp3$`Expression level` <- as.numeric(temp3$`Expression level`) ### 최종 ggplot 테이블
  }
  ####################################
  
  ### ggplot 그래프 그리기 ###
  {
    #초기 그래프 
    p <- ggplot(temp3, aes(y=`Expression level`, x=`Group`, fill=`Group`)) +  
      geom_boxplot(position=position_dodge(1), color = "black", linewidth = 0.5)
    
    # 그래프에 1) 점 표시 / 2) X 축 라벨 지정 / 3) Y 축 범위 및 간격 지정 / 4) 그래프 타이틀 및  X축, Y축 이름 지정 / 
    # 5) 그래프 색지정 / 6) 레전드 제거 및 글씨체 및 글씨크기 설정 / 7), 8) xy축 글씨 크기 및 색깔 위치 지정 / 
    # 9) 축 라인 굵기, 스타일, 색 설정 / 10) 축의 틱 굵기, 길이, 색 설정 / 11) 유의수준 표시 / 12) 그래프 저장
    
    p + 
      geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), 
                   dotsize = 0.5, color = "black", fill = "black") + #1)
      scale_x_discrete(limits=c("HC", "SS")) + #2)
      scale_y_continuous(limits = c(FminY, FmaxY), breaks = seq(FminY, FmaxY, tickgap)) + #3)
      labs(title = i, x="", y = "Expression level\n", face = "italic") + #4)
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
  }
  ############################
  
  ### 그래프 저장 ###
  ggsave(filename = fname1, dpi = 72, width = 150, height = 300, units = 'px') #12)
  
}


############ Correlation Plot ####################
{
  
  library('DGCA')
  
  Sum_Corr6
  
  group <- factor(c(rep('SS', 39), rep('HC', 40)))
  design <- model.matrix(~0+group)
  colnames(design) <- c("SS", "HC")
  str(design)
  
  ddcor_res = ddcorAll(inputMat = NCC2, design = design,
                       compare = c("SS", "HC"))
  
  #BiocManager::install('MEGENA')
  library('MEGENA')
  
  
  
  
  plotCors(NCC2, design, compare = c("HC", "SS"), corrType = "pearson", geneA = "HP", geneB = "TNF",
           oneRow = FALSE, smooth = TRUE, log = FALSE, ylab = NULL,
           xlab = NULL)
  
  SS_Corr22 <- SS_Corr[SS_Corr$corr > 0 & SS_Corr$p.value < 0.05,]
  SS_Corr23 <- SS_Corr22[!is.na(SS_Corr22$corr),]
  dim(SS_Corr23)
  
  
}