### Fig 4d Volcano plot ###

deg <- read.csv(file = "d:/HP_rdata/deg_pbmc.combined4.csv", row.names = 1)




# Volcano plot 그리기
{
  library(ggplot2)
  library(extrafont)
  library('ggrepel')
  
  {
    fname = "MNvsMC.pdf"
    qlf1 <- deg
    IF_genes <- read.csv(file = "d:/GO_term/cytokine_chemokine2.csv")

    qlf1$label <- rownames(qlf1)
    qlf1$IF_label <- NA
    qlf1$IF_label[rownames(qlf1) %in% c(IF_genes$Cytokine, 
                                    IF_genes$Chemokine, "TNF")] <- rownames(qlf1)[qlf1$label %in% c(IF_genes$Cytokine, 
                                                                                                   IF_genes$Chemokine, "TNF")]
    
    qlf1$DE <- "ns"
    qlf1$DE[qlf1$avg_log2FC > 1 & qlf1$p_val_adj < 0.05] <- "UP"
    qlf1$DE[qlf1$avg_log2FC  < -1 & qlf1$p_val_adj < 0.05] <- "DOWN"
    
    #HC vs MC 라서 fc가 + , - 반대로 됨
    qlf1$IF_label[(qlf1$DE == "UP" | qlf1$DE == "ns")] <- NA

    qlf1$log10pval <- -log10(qlf1$p_val_adj)
    
    max(qlf1$log10pval)
    qlf1$log10pval[order(-qlf1$log10pval)]
    
    qlf1[qlf1$log10pval == "Inf" | qlf1$log10pval > 300 ,9] <- 300

    Upgenes <- qlf1[qlf1$DE == "UP", ]
    Downgenes <- qlf1[qlf1$DE == "DOWN", ]
    
    p <- ggplot(data=qlf1, aes(x=-avg_log2FC, y=log10pval, fill=DE, label = IF_label)) +
      geom_point(shape = 21, fill = 'white', stroke = 1, aes(color = DE), size = 2) +
      #geom_text(col = "black", hjust=-.1, size = 3, check_overlap=TRUE) +
      geom_vline(xintercept=c(-1), col="black", linetype = 2, size = 0.5) + 
      geom_vline(xintercept=c(1), col="black", linetype = 2, size = 0.5) +
      geom_hline(yintercept=-log10(0.05), col="black", linetype = 2, size = 0.5) +
      scale_color_manual(values=c("#1E90FF", "gray", "#CD5C5C")) +
      scale_x_continuous(limits=c(0,10), breaks=c(0,seq(1, 10, 3))) +
      scale_y_continuous(limits=c(0,300), breaks=c(seq(0, 300, 100))) +
      theme_classic() +
      theme(plot.title = element_text(family = "Arial Narrow", hjust = 0.5, size = 30, color = "black")) + # theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))   # 글씨체, 글씨 모양, 가운데 정렬, 크기, 색상을 설정합니다.
      theme(legend.position = "none") +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
      theme(text = element_text(size = 16), #6)
            axis.text.x = element_text(size = 16, vjust = -0.5, color = "black"), #7)
            axis.text.y = element_text(size = 16, color = "black"), #8)
            axis.line = element_line(size = 0.5, linetype = "solid", color = "black")) + #9)
      xlab(bquote("log"["2"]~"F.C")) +
      ylab(bquote("-log"["10"]~"p value")) +
      geom_label_repel(size = 3.5, colour="black", max.overlaps = getOption("ggrepel.max.overlaps", default = 30), 
                       box.padding = 1) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    #geom_label_repel(size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
    
    ggsave(filename = fname, dpi = 72, width = 350, height = 300, units = 'px') #12)
  }
}
