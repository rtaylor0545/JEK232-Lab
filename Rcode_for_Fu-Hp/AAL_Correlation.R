### HP correlation ###

table <- read.csv(file = "d:/HP_correlation.csv")
str(table)

table2 <- table[3:42,]
table2$Sample <- table2$Group
table2$Group[1:20] <- "HC"
table2$Group[21:40] <- "SS"

table2_HC <- table2[1:20,]
table2_SS <- table2[21:40,]

temp <- cor.test(table2_SS$AAL, table2_SS$CLEC4E)
EP <- cbind(Rho = temp$estimate, `p value` = temp$p.value)
rownames(EP) <- "AAL_CLEC4E"

if(round(EP[1,2], 4) == 0) {
  pval <- "0.000"
} else {
  pval <- round(EP[1,2], 4)
}

rlabel <- paste0("Rho = ", signif(round(EP[1,1], 4)), '\n', "p = ", pval, " ")

library(ggplot2)
library(extrafont)

ggplot(table2_SS, aes(x = AAL, y = CLEC4E, color= Group)) +
  geom_point(color = "blue", size=2) +
  theme_bw() +
  xlab("AAL") +
  ylab("CLEC4E") +
  geom_smooth(color = "#0066CC", method = "lm") +
  annotate("text", colour="black", size=5, x=Inf, y=Inf, label=rlabel, hjust=1.2, vjust=1.2, family="Arial Narrow") +
  theme(legend.position="none", text = element_text(size = 16, family="Arial Narrow"), #6)
        axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
        axis.text.y = element_text(size = 16, color = "black"), #8)
        axis.line = element_line(size = 0.5, linetype = "solid", color = "black"), #9)
        axis.ticks = element_line(size = 0.5, linewidth = 10, linetype = "solid", color = "black"), 
        plot.title = element_text(face = "italic"))

ggsave(filename = fname2, dpi = 72, width = 300, height = 300, units = 'px') #12)