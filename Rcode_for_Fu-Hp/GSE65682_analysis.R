#### Sepsis Hp Project ####
#### Cohort : GSE154918 <- getGEO로 안됨####
#### Cohort : GSE65682 ####

library(GEOquery)
geo_name <- "GSE65682"

# Seq data preprocessing -------------------------------------------------------
gse <- getGEO(geo_name, 
              GSEMatrix=TRUE, 
              AnnotGPL = T,
              getGPL = TRUE,
              parseCharacteristics = T)
View(gse)
gse <- gse[[1]]
eset <- exprs(gse)
eset[1:3, 1:3] ### RNAseq data confirm 

ginfo <- cbind(gse@featureData@data[["Gene Symbol"]], gse@featureData@data[["ID"]])
ginfo <- as.data.frame(ginfo)

View(gse)
SYMBOL <- gse@featureData@data[["Gene Symbol"]]
View(gse@phenoData@data)

diag <- gse@phenoData@data[["characteristics_ch1.6"]]

names(diag) <- rownames(gse@phenoData@data)

library(stringr)

info <- gse@phenoData@data
head(info)
colnames(info)
info2 <- info[, c("title", "characteristics_ch1.6")]
HC_code <- info2[str_detect(info2$title,'healthy'), ]
SO_code <- info2[str_detect(info2$characteristics_ch1.6,': 0'), ]
SX_code <- info2[str_detect(info2$characteristics_ch1.6,': 1'), ]

dim(HC_code)
dim(SO_code)
dim(SX_code)

rownames(HC_code)
rownames(SO_code)
rownames(SX_code)

eset[,rownames(HC_code)] #HC - 42 person
eset[,rownames(SO_code)] #SO - 365 person
eset[,rownames(HC_code)] #SX - 114 person

fset <- cbind(eset[,rownames(HC_code)], eset[,rownames(SO_code)], eset[,rownames(SX_code)])
head(fset)

write.csv(fset, file = "GSE65682_eset(HC42_SO365_SX114).csv")
write.csv(ginfo, file = "GSE65682_eset(Gene_code).csv")
write.csv(HC_code, file = "GSE65682_eset(HC_code).csv")
write.csv(SO_code, file = "GSE65682_eset(SO_code).csv")
write.csv(SX_code, file = "GSE65682_eset(SX_code).csv")

# Analysis ---------------------------------------------------------------------

library(edgeR)

group <- factor(c(rep(1, 42), rep(2, 365), rep(3, 114)))

y <- DGEList(counts=fset, group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)

rownames(ginfo) <- ginfo$V2
head(ginfo[rownames(qlf.2vs1$table), ])
head(qlf.2vs1$table)

fqlf.2vs1 <- cbind(Gene = ginfo[rownames(qlf.2vs1$table), 1], qlf.2vs1$table)
colnames(fqlf.2vs1)

# To compare 3 vs 1:
qlf.3vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.3vs1)

fqlf.3vs1 <- cbind(Gene = ginfo[rownames(qlf.3vs1$table), 1], qlf.3vs1$table)

# To compare 3 vs 2:
qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))
topTags(qlf.3vs2)

fqlf.3vs2 <- cbind(Gene = ginfo[rownames(qlf.3vs2$table), 1], qlf.3vs2$table)

# The contrast argument in this case requests a statistical test of the null hypothesis that
# coefficient3−coefficient2 is equal to zero.
# To find genes different between any of the three groups:
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)

write.csv(fqlf.2vs1, file = "01_HC_vs_SRO(GSE65682).csv")
write.csv(fqlf.3vs1, file = "02_HC_vs_SRX(GSE65682).csv")
write.csv(fqlf.3vs2, file = "03_SRX_vs_SRO(GSE65682).csv")

######## GO inflammation gene 가져오기 ############
IF_gene <- read.csv(file = "d:/GO_term/Hsa_Inflammation_genes.csv")
length(IF_gene$GOBP_Hsa_Inflammatory.response)
length(IF_gene$HM_Hsa_Inflammatory.response)
temp <- IF_gene$HM_Hsa_Inflammatory.response[!(IF_gene$HM_Hsa_Inflammatory.response == "")]


All_genes <- c(IF_gene$GOBP_Hsa_Inflammatory.response, temp)
length(All_genes)

All_genes <- All_genes[-which(duplicated(All_genes))]
length(All_genes)
All_genes <- c(All_genes, "CLEC4E")


table1 <- fqlf.2vs1[(fqlf.2vs1$Gene %in% All_genes) | fqlf.2vs1$Gene == "CLEC4E",]
table2 <- fqlf.3vs1[(fqlf.3vs1$Gene %in% All_genes) | fqlf.3vs1$Gene == "CLEC4E",]
table3 <- fqlf.3vs2[(fqlf.3vs2$Gene %in% All_genes) | fqlf.3vs2$Gene == "CLEC4E",]

dim(table1)

write.csv(table1, file = "01_SRO_vs_HC(IF_GSE65682).csv")
write.csv(table2, file = "02_SRX_vs_HC(IF_GSE65682).csv")
write.csv(table3, file = "03_SRX_vs_SRO(IF_GSE65682).csv")

########## Correlation with CLEC4E(Mincle) ############

dim(ginfo)
dim(y$counts)
ginfo[ginfo$V1 %in% (All_genes), 2]

head(y$counts)

IFeset <- y$counts[rownames(y$counts) %in% ginfo[ginfo$V1 %in% All_genes,2],]
head(IFeset)
dim(IFeset)

IFeset_HC <- as.data.frame(t(cbind(IFeset[,rownames(HC_code)])))
IFeset_SRO <- as.data.frame(t(cbind(IFeset[,rownames(SO_code)])))
IFeset_SRX <- as.data.frame(t(cbind(IFeset[,rownames(SX_code)])))

dim(IFeset_HC)
dim(IFeset_SRO)
dim(IFeset_SRX)

ginfo[ginfo$V1 == "CLEC4E", 2]

target <- "11745144_a_at" # CLEC4E gene code : 11745144_a_at
target_eset <- IFeset_HC[, target]
count <- dim(IFeset_HC)[2]

#----> IFeset_HC <-----#
for(i in 1:count) {
  temp <- cor.test(target_eset, IFeset_HC[, i])
  EP <- cbind(Rho = temp$estimate, `p value` = temp$p.value)
  rownames(EP) <- colnames(IFeset_HC)[i]
  
  if(i == 1){
    data <- EP 
  } else {
    data <- rbind(data, EP)
  }
  
  if(i == count){
    write.csv(data, file ="Cor_with_CLEC4E_in_HC(GSE65682).csv")
  }
}

data <- cbind(Gene = ginfo[rownames(data), 1], data)
write.csv(data, file ="Cor_with_CLEC4E_in_HC(GSE65682).csv")

target_eset <- IFeset_SRO[, target]

#----> IFeset_SRO <-----#
for(i in 1:count) {
  temp <- cor.test(target_eset, IFeset_SRO[, i])
  EP <- cbind(Rho = temp$estimate, `p value` = temp$p.value)
  rownames(EP) <- colnames(IFeset_SRO)[i]
  
  if(i == 1){
    data <- EP 
  } else {
    data <- rbind(data, EP)
  }
  
  if(i == count){
    write.csv(data, file ="Cor_with_CLEC4E_in_SRO(GSE65682).csv")
  }
}

data <- cbind(Gene = ginfo[rownames(data), 1], data)
write.csv(data, file ="Cor_with_CLEC4E_in_SRO(GSE65682).csv")

target_eset <- IFeset_SRX[, target]

#----> IFeset_SRX <-----#
for(i in 1:count) {
  temp <- cor.test(target_eset, IFeset_SRX[, i])
  EP <- cbind(Rho = temp$estimate, `p value` = temp$p.value)
  rownames(EP) <- colnames(IFeset_SRX)[i]
  
  if(i == 1){
    data <- EP 
  } else {
    data <- rbind(data, EP)
  }
  
  if(i == count){
    write.csv(data, file ="Cor_with_CLEC4E_in_SRX(GSE65682).csv")
  }
}

data <- cbind(Gene = ginfo[rownames(data), 1], data)
write.csv(data, file ="Cor_with_CLEC4E_in_SRX(GSE65682).csv")

##############

#install.packages('ggcorrplot')
library(ggcorrplot)
library(corrplot)
library(dplyr)
library(colorspace)

data <- as.data.frame(data)
data$Rho <- as.numeric(data$Rho)
data$`p value` <- as.numeric(data$`p value`)

data2 <- data[order(-data$Rho), ] # 내림차순 정렬
data3 <- data2[data2$`p value` < 0.05, ] # p value < 0.05 필터
data4 <- data3[-which(duplicated(data3$Gene)), ]

dim(data4)

cor_gene <- c("11745144_a_at", 
  "11726337_a_at", 
  "11745878_x_at", 
  "11754833_a_at", 
  "11719657_a_at", 
  "11728552_a_at", 
  "11725138_at", 
  "11753823_a_at", 
  "11731424_x_at", 
  "11743730_at", 
  "11743636_at", 
  "11729349_at")

data5 <- data4[cor_gene,]

IFeset_HC2.cor <- cor(IFeset_HC[, rownames(data5)])

p.mat <- cor_pmat(IFeset_HC, method = "pearson")
p.mat <- p.mat[rownames(IFeset_HC2.cor),colnames(IFeset_HC2.cor)]

colnames(IFeset_HC2.cor) <- data5$Gene
rownames(IFeset_HC2.cor) <- data5$Gene

#HC
p1 <- ggcorrplot(IFeset_HC2.cor,
                 # type='lower',
                 # hc.order=TRUE,
                 lab=TRUE,
                 outline.color='white',
                 p.mat=p.mat,
                 insig='blank',
                 colors=diverge_hcl(3, palette='Blue Red2'))

#SRO
IFeset_SRO.cor <- cor(IFeset_SRO[, rownames(data5)])

p.mat <- cor_pmat(IFeset_SRO, method = "pearson")
p.mat <- p.mat[rownames(IFeset_SRO.cor),colnames(IFeset_SRX.cor)]

colnames(IFeset_SRO.cor) <- data5$Gene
rownames(IFeset_SRO.cor) <- data5$Gene

p2 <- ggcorrplot(IFeset_SRO.cor,
                 # type='lower',
                 #hc.order=TRUE,
                 lab=TRUE,
                 outline.color='white',
                 p.mat=p.mat,
                 insig='blank',
                 colors=diverge_hcl(3, palette='Blue Red2'))

#SRX
IFeset_SRX.cor <- cor(IFeset_SRX[, rownames(data5)])
p.mat <- cor_pmat(IFeset_SRX, method = "pearson")
p.mat <- p.mat[rownames(IFeset_SRX.cor),colnames(IFeset_SRX.cor)]

colnames(IFeset_SRX.cor) <- data5$Gene
rownames(IFeset_SRX.cor) <- data5$Gene

p3 <- ggcorrplot(IFeset_SRX.cor,
                 # type='lower',
                 #hc.order=TRUE,
                 lab=TRUE,
                 outline.color='white',
                 p.mat=p.mat,
                 insig='blank',
                 colors=diverge_hcl(3, palette='Blue Red2'))

final_P <- p1 + p2 + p3
saveRDS(final_P, file = "Correlation_Prof_Pick_genes(GSE65682).rds")
































########################### old ------------------------------------------------
design <- model.matrix(~diag3)

v = voom(eset, design, plot=TRUE)

# Fit model to data given design
fit = lmFit(v, design)
fit = eBayes(fit)

# Show top genes
topGenes = topTable(fit, coef=ncol(design), number=dim(eset)[1], sort.by="p")

head(topGenes)
topGenes$rn <- rownames(topGenes)

library(dplyr)
topGenes2 <- left_join(topGenes, ginfo, by = c("rn" = "V2"))


# GEO2R -------------------------------------------------------------------

topGenes2 <- readRDS(file = "D:/AML/AML_Rdata/GSE9476_foldchange_table.rds")
str(topGenes2)

# Box plot ----------------------------------------------------------------
final <- eset
colnames(final) <- diag3
final[1:3, 1:3]
final[which(rownames(final) == "203193_at"), colnames(final) == "HC"]
final[which(rownames(final) == "203193_at"), colnames(final) == "AML"]

boxplot(final[which(rownames(final) == "203193_at"), colnames(final) == "HC"], 
        final[which(rownames(final) == "203193_at"), colnames(final) == "AML"])

boxplot(final[which(rownames(final) == "1487_at"), colnames(final) == "HC"], 
        final[which(rownames(final) == "1487_at"), colnames(final) == "AML"])

# overlap -----------------------------------------------------------------

sigfinal <- topGenes2[(topGenes2$adj.P.Val < 0.05), ]
dim(sigfinal)

AML3g <- NULL
AML3g$HR1big <- big_overlap
AML3g$HR1sma <- small_overlap
AML3g$Total <- c(big_overlap, small_overlap)

#saveRDS(AML3g, file = "AML3g_overlap.rds")

AML3g <- readRDS(file = "D:/AML/AML_Rdata/AML3g_overlap.rds")

library(VennDiagram)
foverlap <- calculate.overlap(
  x = list(
    A = AML3g$Total,
    B = sigfinal$SYMBOL
  )
)
View(foverlap)
foverlap$a3

# saveRDS(foverlap$a3, file = "(3g_and_FC)final_overlap_in_AML.rds")
# write.csv(foverlap$a3, file = "(3g_and_FC)final_overlap_in_AML.csv")

toverlap <- readRDS(file = "D:/AML/AML_Rdata/AML3gfc_overlap_FC_table2.rds")

topGenes2[topGenes2$V1 == "SMAD3", ]

topGenes2$label <- ifelse((topGenes2$V1 %in% foverlap$a3), T, F)
topGenes2$label[topGenes2$adj.P.Val > 0.05] <- F
topGenes2$diffexpressed <- ifelse((topGenes2$logFC > 0), "Up-regulated", "Down-regulated")
topGenes2$diffexpressed[topGenes2$adj.P.Val > 0.05] <- NA
topGenes2$marker <- ifelse(topGenes2$label, topGenes2$V1, NA)
head(topGenes2)

topGenes3 <- topGenes2[(topGenes2$V1 %in% foverlap$a3) & !is.na(topGenes2$diffexpressed), ]
topGenes3 <- topGenes3[!duplicated(topGenes3$V1), ]

dim(topGenes3)

#??ġ?? marker ��?? (?? ó��?? ???°͸? ?츮??)
topGenes2$marker[!is.na(topGenes2$marker)][duplicated(topGenes2$marker[!is.na(topGenes2$marker)])] <- NA

#????��?? ?߷?�� (HR ?????ؼ?)
marker <- c("212033_at",
            "204108_at",
            "209999_x_at",
            "221284_s_at",
            "49878_at",
            "209437_s_at",
            "211434_s_at",
            "202476_s_at",
            "220042_x_at",
            "205397_x_at"
)

topGenes2$marker <- NA
topGenes2$marker[topGenes2$rn %in% marker] <- topGenes2$V1[topGenes2$rn %in% marker]

saveRDS(topGenes2, file = "GSE9476_foldchange_table2.rds")
write.csv(topGenes2, file = "GSE9476_foldchange_table2.csv")
saveRDS(topGenes3, file = "AML3gfc_overlap_FC_table2.rds")
write.csv(topGenes3, file = "AML3gfc_overlap_FC_table2.csv")

topGenes2 <- readRDS(file = "D:/AML/AML_Rdata/GSE9476_foldchange_table2.rds")

library('ggrepel')

#aml ????��?? ?ٲٱ? 
ggplot(data=topGenes2, aes(x=-logFC, y=-log10(adj.P.Val), col=diffexpressed, label = marker)) +
  geom_point(aes(alpha= rev(-log10(adj.P.Val)))) + 
  #geom_text(col = "black", hjust=-.1, size = 3, check_overlap=TRUE) +
  #geom_vline(xintercept=c(-1), col="black", linetype = 2, size = 0.5) + 
  #geom_vline(xintercept=c(1), col="black", linetype = 2, size = 0.5) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 0.5) +
  scale_color_manual(values=c("red", "blue", "gray")) +
  theme_classic() +
  theme(plot.title = element_text(family = "serif", hjust = 0.5, size = 30, color = "black")) + # theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))   # ?۾?ü, ?۾? ????, ??? ��??, ũ??, ????�� ??��?մϴ?.
  theme(legend.position = "none") +
  geom_label_repel(size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 30))


# GEO2R data��?? volcano plot ?׸???

table <- read.csv(file = "D:/AML/AML_Rdata/GSE9476_foldchange.csv")
table[1:3, 1:3]
colnames(table)

table$marker <- NA
table$marker[table$logFC > 0 & table$adj.P.Val < 0.05] <- "UP"
table$marker[table$logFC < 0 & table$adj.P.Val < 0.05] <- "DOWN"

View(table)

# Volcano graph ?׸???
ggplot(data=table, aes(x=-logFC, y=-log10(adj.P.Val), col = marker)) +
  geom_point(aes(alpha= rev(-log10(adj.P.Val))), size = 3) + 
  #geom_text(col = "black", hjust=-.1, size = 3, check_overlap=TRUE) +
  #geom_vline(xintercept=c(-1), col="black", linetype = 2, size = 0.5) + 
  #geom_vline(xintercept=c(1), col="black", linetype = 2, size = 0.5) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 0.5) +
  scale_color_manual(values=c("red", "blue", "gray")) +
  theme_classic() +
  theme(plot.title = element_text(family = "serif", hjust = 0.5, size = 30, color = "black")) + # theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))   # ?۾?ü, ?۾? ????, ??? ��??, ũ??, ????�� ??��?մϴ?.
  theme(legend.position = "none") +
  xlab(bquote("log"["2"]~"FoldChange")) +
  ylab(bquote("-log"["10"]~"(adj.P-value)")) +
  theme(
    # LABLES APPEARANCE
    plot.title = element_text(size=14, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=18, face="bold", colour = "black"),    
    axis.title.y = element_text(size=18, face="bold", colour = "black"),    
    axis.text.x = element_text(size=14, face="plain", colour = "black"), 
    # axis.text.y = element_text(size=12,  colour = "black"), # unbold
    axis.text.y = element_text(size=14, face="plain", colour = "black"), # bold
    strip.text.x = element_text(size = 10, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 10, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )
#  geom_label_repel(size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
