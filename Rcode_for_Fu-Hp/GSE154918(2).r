temp <- read.csv(file = "d:/GSE154918_Schughart_Sepsis_200320.txt", sep = " ")
head(temp)

dim(temp)
temp2 <- temp[,10:114]
rownames(temp2) <- temp$ENSEMBL_gene_ID
head(temp2)

library(GEOquery)
geo_name <- "GSE154918"

gse <- getGEO(geo_name, 
              GSEMatrix=TRUE, 
              AnnotGPL = T,
              getGPL = TRUE,
              parseCharacteristics = T)
gse <- gse[[1]]

cohort <- as.data.frame(cbind(ID = gse@phenoData@data[["title"]], Status = gse@phenoData@data[["characteristics_ch1"]]))

library(stringr)
HC_eset <- temp2[, which(str_detect(cohort$Status,'Hlty'))]
SP_eset <- temp2[, which(str_detect(cohort$Status,'Seps_P'))]
SK_eset <- temp2[, which(str_detect(cohort$Status,'Shock_P'))]
SPF_eset <- temp2[, which(str_detect(cohort$Status,'Seps_FU'))]
SKF_eset <- temp2[, which(str_detect(cohort$Status,'Shock_FU'))]

dim(HC_eset)
dim(SP_eset)
dim(SK_eset)
dim(SPF_eset)
dim(SKF_eset)

library(edgeR)

fset <- cbind(HC_eset, SP_eset, SK_eset)
group <- factor(c(rep(1, 40), rep(2, 20), rep(3, 19)))

write.csv(rownames(fset), file = "esbl.csv")
ginfo <- read.csv(file = "esbl.csv")

y <- DGEList(counts=fset, group=group)
# keep <- filterByExpr(y)
# y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

write.csv(fit$counts, file = "GSE154918.csv")
saveRDS(fit$counts, file = "GSE154918_eset.rds")
saveRDS(fit, file = "GSE154918_edgeR_product.rds")

qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)

head(qlf.2vs1$table)
qlf.2vs1$table$Gene <- ginfo$Gene 
qlf.2vs1$table

fqlf.2vs1 <- qlf.2vs1$table
# To compare 3 vs 1:
qlf.3vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.3vs1)

qlf.3vs1$table$Gene <- ginfo$Gene 
qlf.3vs1$table

fqlf.3vs1 <- qlf.3vs1$table
# To compare 3 vs 2:
qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))
topTags(qlf.3vs2)

qlf.3vs2$table$Gene <- ginfo$Gene 
qlf.3vs2$table

fqlf.3vs2 <- qlf.3vs2$table
# The contrast argument in this case requests a statistical test of the null hypothesis that
# coefficient3−coefficient2 is equal to zero.
# To find genes different between any of the three groups:
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)

write.csv(qlf.2vs1$table, file = "01_SRO_vs_HC(GSE154918).csv")
write.csv(qlf.3vs1$table, file = "02_SRX_vs_HC(GSE154918).csv")
write.csv(qlf.3vs2$table, file = "03_SRX_vs_SRO(GSE154918).csv")

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

write.csv(table1, file = "01_SRO_vs_HC(IF_GSE154918).csv")
write.csv(table2, file = "02_SRX_vs_HC(IF_GSE154918).csv")
write.csv(table3, file = "03_SRX_vs_SRO(IF_GSE154918).csv")

########## Correlation with CLEC4E(Mincle) ############

dim(ginfo)
dim(y$counts)
length(All_genes)
ginfo[ginfo$Gene %in% (All_genes), 2]

head(y$counts)

cor_gene2 <- c("MMP25",
               "MMP9",
               "MMP8",
               "CXCL17",
               "TLR8",
               "IL18R1",
               "ALOX5",
               "TNFAIP6",
               "IL17RA",
               "S100A8",
               "S100A9",
               "CLEC4E")

IFeset <- y$counts[rownames(y$counts) %in% ginfo[ginfo$Gene %in% (cor_gene2), 1],]
head(IFeset)
dim(IFeset)

colnames(HC_eset)
colnames(SP_eset)
colnames(SK_eset)

IFeset_HC <- as.data.frame(t(cbind(IFeset[,colnames(HC_eset)])))
IFeset_SRO <- as.data.frame(t(cbind(IFeset[,colnames(SP_eset)])))
IFeset_SRX <- as.data.frame(t(cbind(IFeset[,colnames(SK_eset)])))

dim(IFeset_HC)
dim(IFeset_SRO)
dim(IFeset_SRX)

ginfo[ginfo$Gene == "CLEC4E", 1]

target <- "ENSG00000166523" # CLEC4E gene code : ENSG00000166523
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
    data <- cbind(Gene = ginfo[rownames(data), 2], data)
    write.csv(data, file ="Cor_with_CLEC4E_in_HC(GSE154918).csv")
  }
}


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
    data <- cbind(Gene = ginfo[rownames(data), 2], data)
    write.csv(data, file ="Cor_with_CLEC4E_in_SRO(GSE154918).csv")
  }
}

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
    data <- cbind(Gene = ginfo[rownames(data), 2], data)
    write.csv(data, file ="Cor_with_CLEC4E_in_SRX(GSE154918).csv")
  }
}

##############

#install.packages('ggcorrplot')
library(ggcorrplot)
library(corrplot)
library(dplyr)
library(colorspace)


cor_gene2

OPF <- c(OPF, "CLEC4E")
OPF <- OPF[-1]

data5 <- data4[cor_gene,]

OPF2 <- ginfo[which(ginfo$Gene %in% cor_gene2), ]

IFeset_HC2.cor <- cor(IFeset_HC[, OPF2$ENSEMBL_ID])
p.mat <- cor_pmat(IFeset_HC, method = "pearson")

colnames(IFeset_HC2.cor) <- OPF2[colnames(IFeset_HC2.cor), 2]
rownames(IFeset_HC2.cor) <- colnames(IFeset_HC2.cor)

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
IFeset_SRO.cor <- cor(IFeset_SRO[, OPF2$ENSEMBL_ID])
p.mat <- cor_pmat(IFeset_SRO, method = "pearson")

colnames(IFeset_SRO.cor) <- OPF2[colnames(IFeset_SRO.cor), 2]
rownames(IFeset_SRO.cor) <- colnames(IFeset_SRO.cor)

p2 <- ggcorrplot(IFeset_SRO.cor,
                 # type='lower',
                 #hc.order=TRUE,
                 lab=TRUE,
                 outline.color='white',
                 p.mat=p.mat,
                 insig='blank',
                 colors=diverge_hcl(3, palette='Blue Red2'))

#SRX
IFeset_SRX.cor <- cor(IFeset_SRX[, OPF2$ENSEMBL_ID])
p.mat <- cor_pmat(IFeset_SRX, method = "pearson")

colnames(IFeset_SRX.cor) <- OPF2[colnames(IFeset_SRX.cor), 2]
rownames(IFeset_SRX.cor) <- colnames(IFeset_SRX.cor)

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
