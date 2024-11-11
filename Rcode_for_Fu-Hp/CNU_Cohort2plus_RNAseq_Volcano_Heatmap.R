library(VennDiagram)

HC_overlap <- calculate.overlap(
  x = list(
    HC1 = na.omit(HC1$X),
    HC2 = na.omit(HC2$Gene), 
    HC3 = na.omit(HC3$Gene)
  )
)

View(HC_overlap)

Corr_genes <- HC_overlap[["a5"]]

venn.diagram(
  x = list(
    HC1 = na.omit(HC1$X),
    HC2 = na.omit(HC2$Gene), 
    HC3 = na.omit(HC3$Gene)
  ),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE)

HC1 <- read.csv(file = "d:/corr/HC1.csv")
HC2 <- read.csv(file = "d:/corr/HC2.csv")
HC3 <- read.csv(file = "d:/corr/HC3.csv")

SRO1 <- read.csv(file = "d:/corr/SRO1.csv")
SRO2 <- read.csv(file = "d:/corr/SRO2.csv")
SRO3 <- read.csv(file = "d:/corr/SRO3.csv")

SRX1 <- read.csv(file = "d:/corr/SRX1.csv")
SRX2 <- read.csv(file = "d:/corr/SRX2.csv")
SRX3 <- read.csv(file = "d:/corr/SRX3.csv")

O1 <- read.csv(file = "d:/corr/SRO_vs_HC(1).csv")
O2 <- read.csv(file = "d:/corr/SRO_vs_HC(2).csv")
O3 <- read.csv(file = "d:/corr/SRO_vs_HC(3).csv")

OP1 <- O1[O1$PValue <0.05, ]
OP2 <- O2[O2$PValue <0.05, ]
OP3 <- O3[O3$PValue <0.05, ]

OP_overlap <- calculate.overlap(
  x = list(
    OP1 = na.omit(OP1$X),
    OP2 = na.omit(OP2$Gene), 
    OP3 = na.omit(OP3$Gene)
  )
)

OP_genes <- OP_overlap[["a5"]]
OPF <- Corr_genes[Corr_genes %in% OP_genes]

X1 <- read.csv(file = "d:/corr/SRX_vs_HC(1).csv")
X2 <- read.csv(file = "d:/corr/SRX_vs_HC(2).csv")
X3 <- read.csv(file = "d:/corr/SRX_vs_HC(3).csv")

XP1 <- X1[X1$PValue <0.05, ]
XP2 <- X2[X2$PValue <0.05, ]
XP3 <- X3[X3$PValue <0.05, ]

XP_overlap <- calculate.overlap(
  x = list(
    XP1 = na.omit(XP1$X),
    XP2 = na.omit(XP2$Gene), 
    XP3 = na.omit(XP3$Gene)
  )
)

XP_genes <- XP_overlap[["a5"]]
XPF <- Corr_genes[Corr_genes %in% XP_genes]


head(O1)

################################################################################
library(ggplot2)
library('ggrepel')

x <- as.data.frame(cbind(Gene = O1$X, logFC = O1$logFC, P.Val = O1$PValue, DE = NA))
x2 <- x[x$Gene != "", ]

x2$DE[x2$Gene %in% OPF] <- "Sig"
str(x2)
x2$logFC <- as.numeric(x2$logFC)
x2$P.Val <- as.numeric(x2$P.Val)
x2$DE
x2$label <- NA

count <- which(!is.na(x2$Gene[x2$DE == "Sig"]))

for(i in count) {
  x2$label[i] <- x2$Gene[i]
}

factor(x2$label)

p <- ggplot(data=x2, aes(x=logFC, y=-log10(P.Val), fill=DE, label = label)) +
  geom_point(shape = 21, fill = 'white', stroke = 1, aes(color = DE), size = 2) +
  #geom_text(col = "black", hjust=-.1, size = 4, check_overlap=TRUE) +
  #geom_point(shape = 12, size = 4, stroke = 0.8) +
  #geom_point(shape = 21, size = 3, stroke = 1, aes(alpha= rev(-log10(adj.P.Val)))) +
  #geom_text(col = "black", hjust=-.1, size = 3, check_overlap=TRUE) +
  geom_vline(xintercept=c(-1), col="black", linetype = 2, size = 0.5) + 
  geom_vline(xintercept=c(1), col="black", linetype = 2, size = 0.5) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 0.5) +
  #scale_color_manual(values = c("dodgerblue3", "#FF3362")) +
  #scale_fill_manual(values=c("#FF3362", "#33C8FF", "black")) +
  #scale_color_manual(values=c("black")) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_continuous(limits = c(-10, 10)) +
  theme_classic() + 
  theme(plot.title = element_text(family = "serif", hjust = 0.5, size = 30, color = "black")) + # theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))   # ?۾?ü, ?۾? ????, ??? ��??, ũ??, ????�� ??��?մϴ?.
  theme(legend.position = "none") + 
  geom_label_repel(size = 3.5, colour="black", max.overlaps = getOption("ggrepel.max.overlaps", default = 30), 
                   box.padding = 1) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
#xlab(bquote("log"["2"]~"F.C")) +
#ylab(bquote("-log"["10"]~"P-value")) 


################################################################################

#install.packages('ggcorrplot')
library(ggcorrplot)
library(corrplot)
library(dplyr)
library(colorspace)

cor_gene <- c(OPF, "CLEC4E")

data5 <- HC1[HC1$X %in% OPF,]
rownames(data5) <- data5$X
data5 <- data5[,-1]

IFeset_HC2.cor <- cor(data5[, rownames(data5)])

p.mat <- cor_pmat(IFeset_HC, method = "pearson")
p.mat <- p.mat[rownames(IFeset_HC2.cor),colnames(IFeset_HC2.cor)]

colnames(IFeset_HC2.cor) <- data5$Gene
rownames(IFeset_HC2.cor) <- data5$Gene

#HC
p1 <- ggcorrplot(data5,
                 # type='lower',
                 # hc.order=TRUE,
                 lab=TRUE,
                 outline.color='white',
                 p.mat=data5$p.value,
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

################################################################################
cor_gene2
O1T <- O1[O1$X %in% cor_gene2, ]
O2T <- O2[O2$Gene %in% cor_gene2, ]
O3T <- O3[O3$Gene %in% cor_gene2, ]

write.csv(O1T, file = "O1T.csv")
write.csv(O2T, file = "O2T.csv")
write.csv(O3T, file = "O3T.csv")

X1T <- X1[X1$X %in% cor_gene2, ]
X2T <- X2[X2$Gene %in% cor_gene2, ]
X3T <- X3[X3$Gene %in% cor_gene2, ]

write.csv(X1T, file = "X1T.csv")
write.csv(X2T, file = "X2T.csv")
write.csv(X3T, file = "X3T.csv")

