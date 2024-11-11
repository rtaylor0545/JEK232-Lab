# Fig.6j 

AAL <- read.csv(file = "d:/HP_rdata/qPCR_AAL_IL1B_IL6_table.csv", row.names = 1)
AAL <- AAL[!is.na(AAL$Po.No.),]
dim(AAL)

library(stringr)

str_detect(AAL$Po.No.,'HC')
AAL_H <- AAL[str_detect(AAL$Po.No.,'HC'),]
dim(AAL_H)
AAL_H2 <- AAL_H[!is.na(AAL_H$IL1B),]
dim(AAL_H2)

str_detect(AAL$Po.No.,'gd')
AAL_S <- AAL[str_detect(AAL$Po.No.,'gd'),]
dim(AAL_S)
AAL_S2 <- AAL_S[!is.na(AAL_S$IL1B),]
dim(AAL_S2)

str_detect(AAL$Po.No.,'bd')
AAL_NS <- AAL[str_detect(AAL$Po.No.,'bd'),]
dim(AAL_NS)
AAL_NS2 <- AAL_NS[!is.na(AAL_NS$IL1B),]
dim(AAL_NS2)

# Correlation graph

# HC --------
test <- cor.test(AAL_H2$AAL, AAL_H2$IL1B)
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
AAL_H2$Group <- as.character(AAL_H2$Group)
AAL_H2$Group <- "HC"

library(ggplot2)

H1 <- ggplot(AAL_H2, aes(x = AAL, y = IL1B, color = Group)) +
  geom_point(size=2, alpha = 0.5) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL1B") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('gray','blue','red')) +
  theme(legend.position = "none")

# S --------

test <- cor.test(AAL_S2$AAL, AAL_S2$IL1B)
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
AAL_S2$Group <- as.character(AAL_S2$Group)
AAL_S2$Group <- "HC"

library(ggplot2)

S1 <- ggplot(AAL_S2, aes(x = AAL, y = IL1B, color = Group)) +
  geom_point(size=2, alpha = 0.5) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL1B") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('blue','red'))+
  theme(legend.position = "none")

# NS --------

test <- cor.test(AAL_NS2$AAL, AAL_NS2$IL1B)
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
AAL_NS2$Group <- as.character(AAL_NS2$Group)
AAL_NS2$Group <- "HC"

library(ggplot2)

NS1 <- ggplot(AAL_NS2, aes(x = AAL, y = IL1B, color = Group)) +
  geom_point(size=2, alpha = 0.5) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL1B") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('red'))+
  theme(legend.position = "none")

H1 + S1 + NS1 + guides(fill="none")

######IL6-----------------------------------------------------------------------
# HC --------
test <- cor.test(AAL_H2$AAL, AAL_H2$IL6)
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
AAL_H2$Group <- as.character(AAL_H2$Group)
AAL_H2$Group <- "HC"

library(ggplot2)

H11 <- ggplot(AAL_H2, aes(x = AAL, y = IL6, color = Group)) +
  geom_point(size=2, alpha = 0.5) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL6") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('gray','blue','red')) +
  theme(legend.position = "none")

# S --------

test <- cor.test(AAL_S2$AAL, AAL_S2$IL6)
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
AAL_S2$Group <- as.character(AAL_S2$Group)
AAL_S2$Group <- "HC"

library(ggplot2)

S11 <- ggplot(AAL_S2, aes(x = AAL, y = IL6, color = Group)) +
  geom_point(size=2, alpha = 0.5) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL6") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('blue','red'))+
  theme(legend.position = "none")

# NS --------

test <- cor.test(AAL_NS2$AAL, AAL_NS2$IL6)
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
AAL_NS2$Group <- as.character(AAL_NS2$Group)
AAL_NS2$Group <- "HC"

library(ggplot2)

NS11 <- ggplot(AAL_NS2, aes(x = AAL, y = IL6, color = Group)) +
  geom_point(size=2, alpha = 0.5) +
  theme_bw() +
  xlab("AAL") +
  ylab("IL6") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="red", size=3, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2)  +
  scale_color_manual(values=c('red'))+
  theme(legend.position = "none")

H11 + S11 + NS11 + guides(fill="none")
#--------------------------------