fctable <- read.csv(file = "d:/FC_scRNA.csv", row.names = 1)
pos_fc <- fctable[fctable$log2.Fold.Change. > 0,]
neg_fc <- fctable[fctable$log2.Fold.Change. < 0,]

organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(clusterProfiler)
library(enrichplot)

colnames(pos_fc) <- c("logFC", "Mac", "Mo", "adj.pval", "pval")
colnames(neg_fc) <- c("logFC", "Mac", "Mo", "adj.pval", "pval")

geneList <- ftemp[ftemp$PValue<0.05,]
dim(geneList)

original_gene_list <- neg_fc$logFC
names(original_gene_list) <- rownames(neg_fc)
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
               #nPerm        = 10000,
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
write.csv(final, file = "d:/scRNAseq_KEGG_result(neg).csv")
}