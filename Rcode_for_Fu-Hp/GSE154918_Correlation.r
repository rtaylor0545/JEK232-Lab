### UBXD1 껴서 보기, 두번째는 Autophagy gene 껴서 보기. 

saveRDS(ginfo, file = "ENSEMBL_ID.rds")
saveRDS(HC_eset, file = "HC_eset.rds") # Healthy control 
saveRDS(SK_eset, file = "SK_eset.rds") # Sepsis_shock
saveRDS(SP_eset, file = "SP_eset.rds") # Sepsis

ginfo <- readRDS(file = "ENSEMBL_ID.rds")
HC_eset <- readRDS(file = "HC_eset.rds") # Healthy control 
SK_eset <- readRDS(file = "SK_eset.rds") # Sepsis_shock
SP_eset <- readRDS(file = "SP_eset.rds") # Sepsis

head(HC_eset)
cor_gene3 <- c(cor_gene2, "UBXN6")
gene_table <- ginfo[which(ginfo$Gene %in% cor_gene3), ]

cHC_eset <- as.data.frame(t(HC_eset[gene_table$ENSEMBL_ID,]))
cSP_eset <- as.data.frame(t(SP_eset[gene_table$ENSEMBL_ID,]))
cSK_eset <- as.data.frame(t(SK_eset[gene_table$ENSEMBL_ID,]))

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
HC_p.mat <- corr.p(temp$r, n=13)
HC_p.mat$p
colnames(rHC.cor) <- gene_table[colnames(rHC.cor), 2]
rownames(rHC.cor) <- colnames(rHC.cor)
colnames(HC_p.mat$p) <- colnames(rHC.cor)
rownames(HC_p.mat$p) <- colnames(rHC.cor)

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
SP_p.mat <- corr.p(rSP.cor, n=13)
SP_p.mat$p
colnames(rSP.cor) <- gene_table[colnames(rSP.cor), 2]
rownames(rSP.cor) <- colnames(rSP.cor)
colnames(SP_p.mat$p) <- colnames(rSP.cor)
rownames(SP_p.mat$p) <- colnames(rSP.cor)

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
SK_p.mat <- corr.p(temp$r, n=13)
SK_p.mat$p
colnames(rSK.cor) <- gene_table[colnames(rSK.cor), 2]
rownames(rSK.cor) <- colnames(rSK.cor)
colnames(SK_p.mat$p) <- colnames(rSK.cor)
rownames(SK_p.mat$p) <- colnames(rSK.cor)

ph_SK <- ggcorrplot(rSK.cor,
                    #               type='lower',
                    #SK.order=TRUE,
                    lab=TRUE,
                    outline.color='white',
                    #               p.mat=SK_p.mat$p,
                    #              insig='blank',
                    colors=diverge_hcl(3, palette='Blue Red2'))


final_P <- ph_HC + ph_SP + ph_SK