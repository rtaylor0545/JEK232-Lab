### Fig2 venn diagram ###
### 01. Preprocessing ###
coh0 <- read.csv(file = "d:/HP_rdata/NewCNU_SSvsHC.csv", row.names = 1) #New CNUH
coh1 <- read.csv(file = "d:/HP_rdata/GSE232753_SSvsHC.csv", row.names = 1) #GSE232753
coh2 <- read.csv(file = "d:/HP_rdata/GSE154918_FC_SSvsHC.csv", row.names = 1) #GSE154918

genes <- read.csv(file = "d:/GO_term/Hsa_Inflammation_genes.csv")
genes <- colnames(read.csv(file = "d:/GO_term/innate_immune_response2.csv"))
genes <- read.csv(file = "d:/GO_term/KEGG_IL17_TNF_JAKSTAT_CLR_HIF1_NFkB_NLR.csv")
genes00 <- c(genes$IL17, genes$TNF)
genes01 <- genes00[!duplicated(genes00)]
genes02 <- genes01[!is.na(genes01)]
genes03 <- genes02[order(genes02)]
genes2 <- c(genes$HM_Hsa_Inflammatory.response, genes$GOBP_Hsa_Inflammatory.response)
genes3 <- genes2[!duplicated(genes2)]
genes4 <- genes3[genes3 != ""]
length(genes4)
genes5 <- c(genes$HM_Hsa_Inflammatory.response)
genes6 <- genes5[genes5 != ""]
genes7 <- genes$GOBP_Hsa_Inflammatory.response
genes8 <- genes7[genes7 != ""]
cgenes0 <- read.csv(file = "d:/GO_term/cytokine_chemokine2.csv")

write.csv(genes4, file = "d:/GO_term/HM_GOBP_Hsa_inflammation_genes.csv")
write.csv(genes6, file = "d:/GO_term/HM_Hsa_inflammation_genes.csv")

coh01 <- coh0[order(-coh0$logFC),]
dim(coh01)

coh02 <- coh01[rownames(coh01) %in% genes03,]
dim(coh02)

coh03 <- coh02[coh02$logFC > 1,]
dim(coh03)

library(extrafont)
library("pheatmap")
library("RColorBrewer")

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Blues")[7]))(19))

dev.off()

eset <- read.csv("C:/Users/User/Desktop/01_Sepsis_Hp_Project/NewRNAseq/SSvsHC_MyAnalysis/log2Eset.csv", row.names = 1)
eset2 <- eset[,c(25:36,1:24)]
heatmap_eset <- eset[rownames(coh03),]
heatmap_eset2 <- cbind(heatmap_eset[,25:36], heatmap_eset[,1:24])
dim(heatmap_eset2)
which(rownames(heatmap_eset2) %in% c("MMP13", "CX3CL1", "S100A7", "VEGFD"))
heatmap_eset3 <- heatmap_eset2[!(1:40 %in% which(rownames(heatmap_eset2) %in% c("MMP13", "CX3CL1", "S100A7", "VEGFD"))),]
dim(heatmap_eset3)
FUT_eset <- eset[FUT_target,]
FUT_eset2 <- cbind(FUT_eset[,25:36], FUT_eset[,1:24])
FUT_eset3 <- FUT_eset2[order(rownames(FUT_eset2)),]
FUT_target <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9", "FUT10", "FUT11", "POFUT1", "POFUT2")
Group <- rep(c(rep("HC", 12), rep("SS", 24)), 10)

for(i in 1:length(FUT_target)) {
  if(i == 1) {
    temp <- unlist(as.vector(FUT_eset3[i,]))
  } else {
    temp <- c(temp, unlist(as.vector(FUT_eset3[i,])))
  }
}
temp  
length(FUT_target)
length(temp)

temp2 <- ggFUT[1:36,]
temp2$target2 <- "HP"
temp2$temp <-  unlist(as.vector(eset2["HP",]))
final <- rbind(ggFUT, temp2)

ggFUT <- as.data.frame(cbind(Group, target2, temp))
ggFUT$temp <- as.numeric(ggFUT$temp)
ggFUT$target2 <- factor(ggFUT$target2, levels = FUT_target)

install.packages("ggthemes")

p <- ggplot(data=ggFUT, aes(x=target2, y=temp, fill = Group)) + 
  geom_boxplot(width = 10)+ ylim(c(2.5,20)) +
  ggthemes::theme_base()  
  geom_point(data=ggFUT, aes(y=temp, x=target2, fill = Group)) +
  theme_classic() + ylim(c(0,25))

target <- rep(FUT_target, 36)
target2 <- target[order(target)]

display.brewer.all()

cols <-  rev(colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[7], 
                                brewer.pal(n = 9, name = "Greys")[3], 
                                brewer.pal(n = 9, name = "Blues")[7]))(19))

cols1 <- rev(colorRampPalette(c("#F40303", "#FFFFFF"))(9))
cols2 <- colorRampPalette(c("#0006FF", "#FFFFFF"))(9)

cols <- c(cols2, cols1)
cols3 <- cols[-10]
cols <- colorRampPalette(cols3)(100)

rownames(heatmap_eset3)
except <- c(12, 13, 14, 15, 18, 21, 25, 32, 35) 
al <- c(1:36)
!(al %in% except)
heatmap_eset4 <- heatmap_eset3[!(al %in% except),]

a <- pheatmap(heatmap_eset4, 
              scale = "row", 
              legend = T, 
              color = cols,
              border_color = F, 
              gaps_col = c(12),
              #clustering_distance_rows = "euclidean",
              fontsize = 7,
              #gaps_row = 1:dim(mf7)[1],
              cellwidth = 7.5, 
              cellheight = 10,
              legend_breaks = c(-1, 0, 1),
              breaks = seq(-1, 1, 2/100),
              show_colnames = F,
              show_rownames = T,
              cluster_rows = F,
              cluster_cols = F, 
              angle_col = 0, 
              family = "Arial Narrow")


################################ Correlation ###################################
test <- cor.test(final$temp[final$target2 == "HP"], final$temp[final$target2 == "FUT3"])

HP <- final$temp[final$target2 == "HP"]
FUT3 <- final$temp[final$target2 == "FUT3"]
FUT4 <- final$temp[final$target2 == "FUT4"]
FUT5 <- final$temp[final$target2 == "FUT5"]
FUT6 <- final$temp[final$target2 == "FUT6"]
FUT7 <- final$temp[final$target2 == "FUT7"]
FUT9 <- final$temp[final$target2 == "FUT9"]
FUT10 <- final$temp[final$target2 == "FUT10"]
FUT11 <- final$temp[final$target2 == "FUT11"]
POFUT1 <- final$temp[final$target2 == "POFUT1"]
POFUT2 <- final$temp[final$target2 == "POFUT2"]

test <- cor.test(HP, POFUT2)
test$p.value
test$estimate
EP <- cbind(Rho = test$estimate, `p value` = test$p.value)
rownames(EP) <- "IL1B"
str(test)

if(round(EP["IL1B",2], 4) == 0) {
  pval <- "0.000"
} else {
  pval <- round(test$p.value, 4)
}

rlabel <- paste0("Rho = ", signif(round(EP[1], 3)), '\n', "p = ", pval, " ")

str(final)
#table2$Group <- as.character(table2$Group)

#final2 <- final[final$target2 %in% c("POFUT2", "HP"),]
final3 <- as.data.frame(cbind(final$Group[1:36], HP, POFUT2))
colnames(final3)[1] <- "Group"
final3$HP <- as.numeric(final3$HP)
final3$POFUT2 <- as.numeric(final3$POFUT2)

ppPOFUT2 <- ggplot(final3, aes(x = HP, y = POFUT2, color = Group)) +
  geom_point(size=2, alpha = 0.7) +
  theme_bw() +
  xlab("HP") +
  ylab("POFUT2") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('gray','red'))

ppFUT3 + theme(legend.position = "none") + ppFUT4 + theme(legend.position = "none")+ 
  ppFUT5 + theme(legend.position = "none") + ppFUT6 + theme(legend.position = "none") + 
  ppFUT7 + theme(legend.position = "none") + ppFUT9 + theme(legend.position = "none") + 
  ppFUT10  + theme(legend.position = "none")+ ppFUT11  + theme(legend.position = "none")+ 
  ppPOFUT1 + theme(legend.position = "none") + ppPOFUT2 + theme(legend.position = "none")

AALcor <- read.csv(file = "d:/HP_rdata/HPtype_AAL_Corr.csv", row.names = 1)

test <- cor.test(AALcor$AAL, AALcor$IL1B)
test <- cor.test(AALcor$AAL, AALcor$IL6)
test <- cor.test(AALcor$AAL, AALcor$TNF)
test$p.value
test$estimate
EP <- cbind(Rho = test$estimate, `p value` = test$p.value)
rownames(EP) <- "IL6"
str(test)

if(round(EP["IL6",2], 4) == 0) {
  pval <- "0.000"
} else {
  pval <- round(test$p.value, 4)
}

rlabel <- paste0("Rho = ", signif(round(EP[1], 3)), '\n', "p = ", pval, " ")
#table2$Group <- as.character(table2$Group)

str(AALcor)

cIL1B <- ggplot(AALcor, aes(x = AAL, y = IL1B, color = Group)) +
  geom_point(size=2, alpha = 0.7) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL1B") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('gray','red'))+ theme(legend.position = "none")

cIL6 <- ggplot(AALcor, aes(x = AAL, y = IL6, color = Group)) +
  geom_point(size=2, alpha = 0.7) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL6") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('gray','red'))+ theme(legend.position = "none")

cTNF <- ggplot(AALcor, aes(x = AAL, y = TNF, color = Group)) +
  geom_point(size=2, alpha = 0.7) +
  theme_bw() +
  xlab("AAL") +
  ylab("TNF") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('gray','red')) + theme(legend.position = "none")

cIL1B + cIL6 + cTNF

library(gridExtra)
grid.arrange(cIL1B, cIL6, cTNF, ncol=1)

#########################################################################################

coh11 <- coh1[coh1$logFC > 0,]
rownames(coh11) <- coh11$Gene.Symbol
dim(coh11)

coh21 <- coh2[coh2$logFC > 0,]
rownames(coh21) <- coh21$Gene
dim(coh21)

fu12 <- c(rownames(coh11), coh21$Gene)
fu13 <- fu12[!duplicated(fu12)]


coh00 <- coh0[rownames(coh0) %in% genes02,]
dim(coh00)

coh11 <- coh1[coh1$Gene.Symbol %in% genes02,]
dim(coh11)

coh1[coh1$Gene.Symbol == "IL6",]

coh22 <- coh2[coh2$Gene %in% genes02,]
dim(coh22)

coh000 <- coh00[coh00$logFC > 0,]
dim(coh000)

coh111 <- coh11[coh11$logFC > 0,]
rownames(coh111) <- coh111$Gene.Symbol
dim(coh111)

coh222 <- coh22[coh22$logFC > 0,]
rownames(coh222) <- coh222$Gene
dim(coh222)

### 02. Venn diagram ###
# Venn diagram 그리기
#---Venn diagram---#
{
    library(VennDiagram)
    library(RColorBrewer)
    
    x = list(
      `0` = rownames(coh01),
      `1` = rownames(coh11), 
      `2` = genes8
    )
    
    Up_overlap <- calculate.overlap(x)
    
    write.csv(Up_overlap$a5, file = "Up_overlapX.csv")
    
    y = list(
      `0` = rownames(coh01),
      `1` = coh21$Gene[!is.na(coh21$Gene)], 
      `2` = genes4
    )
    
    Up_overlap2 <- calculate.overlap(y)
    write.csv(Up_overlap2$a5, file = "Up_overlapY.csv")
    
    library(RColorBrewer)
    coln <- length(x)
    venn.diagram(x,
                 filename = "fig2_VD_plot(3).tiff",   # 저장할 파일명을 설정합니다.
                 col = "white",   # 벤 다이어그램 테두리 색상을 설정합니다.
                 lty = 1,   # 벤 다이어그램 테두리 선 모양을 설정합니다.
                 fill = c(brewer.pal(n = 9, name = 'Set1')[1:3]),   # 각 벤 다이어그램 내부 색상을 설정합니다.
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
    
    venn.diagram(x,
                 filename = "fig2_VD_plot(4).tiff",   # 저장할 파일명을 설정합니다.
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
    venn.diagram(y,
                 filename = "fig2_VD_plot(5).tiff",   # 저장할 파일명을 설정합니다.
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
    
    venn.diagram(y,
                 filename = "fig2_VD_plot(6).tiff",   # 저장할 파일명을 설정합니다.
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
}
  