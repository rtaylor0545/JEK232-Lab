### UBXD1 껴서 보기, 두번째는 Autophagy gene 껴서 보기. 

eset <- read.csv(file = "eset.csv")
rownames(eset) <- eset$X
eset <- eset[,-1]
saveRDS(eset, file = "GSE65682_eset.rds")

ginfo <- read.csv(file ="GSE65682_eset(Gene_code).csv")
rownames(ginfo) <- ginfo$V2
ginfo <- ginfo[,-1]
ginfo2 <- ginfo[ginfo$V1 %in% cor_gene3, ]

ginfo3 <- ginfo2[!duplicated(ginfo2$V1), ]

#(HC42_SO365_SX114)
HC_eset <- eset[,1:42]
SP_eset <- eset[,43:407]
SK_eset <- eset[,408:521]

saveRDS(ginfo3, file = "Picked_Gene_annotation.rds") # Picked Gene ID and Symbol
saveRDS(ginfo, file = "Gene_annotation.rds") # Gene ID and Symbol
saveRDS(HC_eset, file = "HC_eset.rds") # Healthy control 
saveRDS(SP_eset, file = "SP_eset.rds") # Sepsis
saveRDS(SK_eset, file = "SK_eset.rds") # Sepsis_shock

HC_eset <- readRDS(file = "HC_eset.rds") # Healthy control 
SP_eset <- readRDS(file = "SP_eset.rds") # Sepsis
SK_eset <- readRDS(file = "SK_eset.rds") # Sepsis_shock

cor_gene3 # 원하는 유전자 삽입

cHC_eset <- as.data.frame(t(HC_eset[ginfo3$V2,]))
cSP_eset <- as.data.frame(t(SP_eset[ginfo3$V2,]))
cSK_eset <- as.data.frame(t(SK_eset[ginfo3$V2,]))

#------------------------------------------------------------------------------#

#install.packages("psych")
library(psych)
library(ggcorrplot)
library(corrplot)
library(dplyr)
library(colorspace)

#HC
HC.cor <- corr.test(cHC_eset)
rHC.cor <- HC.cor$r

colnames(rHC.cor) <- ginfo3[ginfo3$V2,1]
rownames(rHC.cor) <- colnames(rHC.cor)

ph_HC <- ggcorrplot(rHC.cor,
                    #               type='lower',
                    #hc.order=TRUE,
                    lab=TRUE,
                    outline.color='white',
                    #               p.mat=HC_p.mat$p,
                    #              insig='blank',
                    colors=diverge_hcl(3, palette='Blue Red2'))

#SP (Seps_P)
SP.cor <- corr.test(cSP_eset)
rSP.cor <- SP.cor$r

colnames(rSP.cor) <- ginfo3[ginfo3$V2,1]
rownames(rSP.cor) <- colnames(rSP.cor)

ph_SP <- ggcorrplot(rSP.cor,
                    # type='lower',
                    #hc.order=TRUE,
                    lab=TRUE,
                    outline.color='white',
                    #            p.mat=SP_p.mat$p,
                    #           insig='blank',
                    colors=diverge_hcl(3, palette='Blue Red2'))

#SK
SK.cor <- corr.test(cSK_eset)
rSK.cor <- SK.cor$r

colnames(rSK.cor) <- ginfo3[ginfo3$V2,1]
rownames(rSK.cor) <- colnames(rSK.cor)

ph_SK <- ggcorrplot(rSK.cor,
                    #               type='lower',
                    #SK.order=TRUE,
                    lab=TRUE,
                    outline.color='white',
                    #               p.mat=SK_p.mat$p,
                    #              insig='blank',
                    colors=diverge_hcl(3, palette='Blue Red2'))


final_P <- ph_HC + ph_SP + ph_SK
