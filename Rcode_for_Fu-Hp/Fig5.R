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
g656_3 <- g656_2$gene[!duplicated(g656_2$gene)]
length(g656_3)
g154_2 <- g154[(g154$logFC >0 & g154$PValue < 0.05),]
g154_3 <- g154_2$Gene[!is.na(g154_2$Gene)]
length(g154_3)
scRNA2 <- scRNA[(scRNA$log2FC > 0 & scRNA$p_val < 0.05),]
scRNA3 <- scRNA2$Gene[!duplicated(scRNA2$Gene)]
length(scRNA3)

CNU_CL <- CNU_fc["CLEC4E",]
G232_CL <- g232[g232$Gene.Symbol %in% "CLEC4E",]
G636_CL <- g656[g656$gene %in% "CLEC4E",]
G154_CL <- g154[g154$Gene %in% "CLEC4E",]

x <- list(
  `C` = CNU,
  `G2` = g232_3,
  `G6` = g656_3,
  `GOiir` = iir, 
  `sc` = scRNA3
)

library(VennDiagram)
library(RColorBrewer)

coln <- length(x)

venn.diagram(x,
             filename = "fig5_VD_plot(5).tiff",   # 저장할 파일명을 설정합니다.
             col = "white",   # 벤 다이어그램 테두리 색상을 설정합니다.
             lty = 1,   # 벤 다이어그램 테두리 선 모양을 설정합니다.
             fill = c(brewer.pal(n = 9, name = 'Set1')[1:coln]),   # 각 벤 다이어그램 내부 색상을 설정합니다.
             alpha = 0.20,   # 벤 다이어그램의 알파 투명도를 설정합니다.
             cat.col = "black",   # 각 벤 다이어그램의 명칭에 대한 색상을 설정합니다.
             cat.cex = 1,   # 각 벤 다이어그램의 명칭에 대한 글자 크기를 설정합니다.
             cat.fontface = "bold",   # 각 벤 다이어그램의 명칭에 대한 글자를 볼드체로 설정합니다.
             margin = 0.1,   #  벤 다이어그램 주위의 공간을 설정합니다.
             scaled = F, 
             sub.cex = 1, 
             main.cex = 1, 
             cex = 1
)

overlap <- calculate.overlap(x)
View(overlap)

x2 <- list(
  `C` = CNU,
  `G2` = g232_3,
  `G6` = g656_3,
  `sc` = scRNA3
)

library(VennDiagram)
library(RColorBrewer)

coln <- length(x2)

venn.diagram(x2,
             filename = "fig5_VD_plot(4).tiff",   # 저장할 파일명을 설정합니다.
             col = "white",   # 벤 다이어그램 테두리 색상을 설정합니다.
             lty = 1,   # 벤 다이어그램 테두리 선 모양을 설정합니다.
             fill = c(brewer.pal(n = 9, name = 'Set1')[1:coln]),   # 각 벤 다이어그램 내부 색상을 설정합니다.
             alpha = 0.20,   # 벤 다이어그램의 알파 투명도를 설정합니다.
             cat.col = "black",   # 각 벤 다이어그램의 명칭에 대한 색상을 설정합니다.
             cat.cex = 1,   # 각 벤 다이어그램의 명칭에 대한 글자 크기를 설정합니다.
             cat.fontface = "bold",   # 각 벤 다이어그램의 명칭에 대한 글자를 볼드체로 설정합니다.
             margin = 0.1,   #  벤 다이어그램 주위의 공간을 설정합니다.
             scaled = F, 
             sub.cex = 1, 
             main.cex = 1, 
             cex = 1
)

overlap2 <- calculate.overlap(x2)
View(overlap2)
a18 <- overlap$a18
length(overlap)

library(plyr)
data.list.df <- ldply(overlap, function(x) paste(x, collapse = " "))
View(data.list.df)
write.csv(data.list.df, file = "venn_diagram.csv")

View(final)

View(overlap)
str(overlap)
write.lis(overlap, file = "d:/venn_diagram.csv")
