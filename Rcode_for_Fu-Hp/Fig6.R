### Venn diagram ###
CNU_fc <- read.csv("C:/Users/User/Desktop/01_Sepsis_Hp_Project/NewRNAseq/SSvsHC_MyAnalysis/SSvsHC.csv", row.names = 1)
g232 <- read.csv("d:/HP_rdata/GSE232753_SSvsHC.csv", row.names = 1)
g656 <- read.csv("C:/Users/User/Desktop/01_Sepsis_Hp_Project/GSE65682/231228/GSE65682_FC_table.csv", row.names = 1)
g154 <- read.csv("C:/Users/User/Desktop/01_Sepsis_Hp_Project/GSE154918/GSE154918_FC_SSvsHC.csv", row.names = 1)
iir <- read.csv("d:/GO_term/innate_immune_response2.csv")
iir <- colnames(iir)
scRNA <- read.csv("d:/HP_rdata/scRNA_deg(MAvsMO).csv")


CNU <- rownames(CNU_fc[(CNU_fc$logFC > 0 & CNU_fc$PValue < 0.05),])
g232_2 <- g232[(g232$logFC > 0 & g232$PValue < 0.05),]
dim(g232_2)
g232_3 <- g232_2$Gene.Symbol
g656_2 <- g656[g656$DE == "UP" & g656$PValue < 0.05,]
dim(g656_2)
g154_2 <- g154[(g154$logFC >0 & g154$PValue < 0.05),]
dim(g154_2)
scRNA2 <- scRNA[(scRNA$log2FC > 0 & scRNA$p_val < 0.05),]
dim(scRNA2)
scRNA3 <- scRNA2$Gene

CNU_CL <- CNU_fc["CLEC4E",]
G232_CL <- g232[g232$Gene.Symbol %in% "CLEC4E",]
G636_CL <- g656[g656$gene %in% "CLEC4E",]
G154_CL <- g154[g154$Gene %in% "CLEC4E",]

x <- list(
  `C` = CNU,
  `G2` = g232_3,
  `G6` = g656_2,
  `GOiir` = iir, 
  `sc` = scRNA3
)


library(gplots)
v.table <- venn(x)
v.table <- venn(x)

library(RColorBrewer)

  coln <- length(x)
  
  venn.diagram(x,
                    filename = "fig6_VD_plot.tiff",   # 저장할 파일명을 설정합니다.
                    col = "white",   # 벤 다이어그램 테두리 색상을 설정합니다.
                    lty = 1,   # 벤 다이어그램 테두리 선 모양을 설정합니다.
                    fill = c(brewer.pal(n = 9, name = 'Set1')[1:coln]),   # 각 벤 다이어그램 내부 색상을 설정합니다.
                    alpha = 0.20,   # 벤 다이어그램의 알파 투명도를 설정합니다.
                    cat.col = "black",   # 각 벤 다이어그램의 명칭에 대한 색상을 설정합니다.
                    cat.cex = 0,   # 각 벤 다이어그램의 명칭에 대한 글자 크기를 설정합니다.
                    cat.fontface = "bold",   # 각 벤 다이어그램의 명칭에 대한 글자를 볼드체로 설정합니다.
                    margin = 0.1,   #  벤 다이어그램 주위의 공간을 설정합니다.
                    scaled = F, 
               sub.cex = 0, 
               main.cex = 0, 
               cex = 0
  )
  
  overlap <- calculate.overlap(x)
  saveRDS(overlap, file = "Overlap_List.rds")
