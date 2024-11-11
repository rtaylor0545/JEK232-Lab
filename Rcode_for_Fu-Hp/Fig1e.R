### Fig1e ###

NCC <- read.csv(file = "d:/NewCNU_Cohort.csv", row.names = 1)
NCC[1:3, 1:3]
SSset <- NCC[,1:24]
HCset <- NCC[,25:36]

SSset[1:3, 1:3]
HCset[1:3, 1:3]

tSSset <- as.data.frame(t(SSset))
tHCset <- as.data.frame(t(HCset))

tSSset[1:3, 1:3]
tHCset[1:3, 1:3]



### Corr 관련 피겨 만들기 ###
{
  
  Basetarget <- "HP"
  
  # HC 그룹
  {
    for(i in 1:dim(tHCset)[2]) {
      if(i == 1) {
        temp <- cor.test(tHCset[,Basetarget], tHCset[,i])
        ctable <- as.data.frame(cbind(corr = temp$estimate, p.value = temp$p.value))
        rownames(ctable) <- colnames(tHCset)[i]
        
      } else {
        temp <- cor.test(tHCset[,Basetarget], tHCset[,i])
        ttable <- as.data.frame(cbind(corr = temp$estimate, p.value = temp$p.value))
        ctable <- as.data.frame(rbind(ctable, ttable))
        rownames(ctable)[i] <- colnames(tHCset)[i]
      }
    }
    
    write.csv(ctable, file = "HC_HP_Corr.csv")
    saveRDS(ctable, file = "HC_HP_Corr.rds")
  }  
  
  
  # SS 그룹
  {
    for(i in 1:dim(tSSset)[2]) {
      if(i == 1) {
        temp <- cor.test(tSSset[,Basetarget], tSSset[,i])
        ctable <- as.data.frame(cbind(corr = temp$estimate, p.value = temp$p.value))
        rownames(ctable) <- colnames(tSSset)[i]
        
      } else {
        temp <- cor.test(tSSset[,Basetarget], tSSset[,i])
        ttable <- as.data.frame(cbind(corr = temp$estimate, p.value = temp$p.value))
        ctable <- as.data.frame(rbind(ctable, ttable))
        rownames(ctable)[i] <- colnames(tSSset)[i]
      }
    }
    
    write.csv(ctable, file = "SS_HP_Corr.csv")
    saveRDS(ctable, file = "SS_HP_Corr.rds")
  }  
  
  # 유전자만 추리기
  {
    HC_Corr <- readRDS(file = "HC_HP_Corr.rds")
    SS_Corr <- readRDS(file = "SS_HP_Corr.rds")
    
    Sum_Corr <- cbind(HC_Corr, SS_Corr)
    colnames(Sum_Corr) <- c("HC_Corr", "HC_pval", "SS_Corr", "SS_pval")
    
    Sum_Corr2 <- Sum_Corr[!is.na(Sum_Corr$HC_Corr | Sum_Corr$SS_Corr),]
    
    head(Sum_Corr2)
    
    genes <- read.csv(file = "d:/GO_term/IF_TNF.csv")
    
    Sum_Corr3 <- Sum_Corr2[rownames(Sum_Corr2) %in% c(genes$IF_Gene, genes$TNFA_Gene),]
    
    dim(Sum_Corr2)
    saveRDS(Sum_Corr2, file = "Whole_genes_Correlation.rds")
    saveRDS(Sum_Corr3, file = "IR_TNF_genes_Correlation.rds")
    write.csv(Sum_Corr3, file = "IR_TNF_genes_Correlation.csv")
    
    Sum_Corr4 <- Sum_Corr3[(Sum_Corr3$HC_pval > 0.05 & Sum_Corr3$SS_pval < 0.05),]
    Sum_Corr5 <- Sum_Corr4[Sum_Corr4$SS_Corr > 0,]
    
    saveRDS(Sum_Corr5, file = "IR_TNF_genes_Correlation(onlySS_pval).rds")
    write.csv(Sum_Corr3, file = "IR_TNF_genes_Correlation.csv")
    
    
    
    rownames(Sum_Corr2)
    
    Sum_Corr6 <- Sum_Corr3[(Sum_Corr3$SS_Corr > 0 & Sum_Corr3$SS_pval < 0.05),]
    
    dim(Sum_Corr6)
    
    Sum_Corr["CLEC4E",]
    

  }

}

Sum_Corr6 <- Sum_Corr3[(Sum_Corr3$SS_Corr > 0 & Sum_Corr3$SS_pval < 0.25),]
rownames(Sum_Corr6)
genes <- read.csv(file = "d:/GO_term/IF_TNF.csv")
fc <- read.csv(file = "d:/SSvsHC.csv", row.names = 1)
fc[1:3, 1:4]


###
#BiocManager::install('DGCA')
library('DGCA')

Sum_Corr6

group <- factor(c(rep('SS', 24), rep('HC', 12)))
design <- model.matrix(~0+group)
colnames(design) <- c("HC", "SS")
str(design)

ddcor_res = ddcorAll(inputMat = NCC2, design = design,
                     compare = c("SS", "HC"))

#BiocManager::install('MEGENA')
library('MEGENA')




plotCors(NCC2, design, compare = c("HC", "SS"), corrType = "pearson", geneA = "HP", geneB = "TNF",
         oneRow = FALSE, smooth = TRUE, log = FALSE, ylab = NULL,
         xlab = NULL)

SS_Corr22 <- SS_Corr[SS_Corr$corr > 0 & SS_Corr$p.value < 0.05,]
SS_Corr23 <- SS_Corr22[!is.na(SS_Corr22$corr),]
dim(SS_Corr23)
####### Correlation Net Plot #######

genes2 <- read.csv(file = "d:/GO_term/Hsa_Inflammation_genes.csv")
genes2$GOBP_Hsa_Inflammatory.response

{
  library(MEGENA)
  
  # input parameters
  n.cores <- 2; # number of cores/threads to call for PCP
  doPar <-TRUE; # do we want to parallelize?
  method = "pearson" # method for correlation. either pearson or spearman. 
  FDR.cutoff = 0.25 # FDR threshold to define significant correlations upon shuffling samples. 
  module.pval = 0.25 # module significance p-value. Recommended is 0.05. 
  hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
  cor.perm = 100; # number of permutations for calculating FDRs for all correlation pairs. 
  hub.perm = 1000; # number of permutations for calculating connectivity significance p-value. 
  
  # annotation to be done on the downstream
  annot.table=NULL
  id.col = 1
  symbol.col= 2
  ###########
  NCC3 <- SSset[rownames(SSset) %in% c(genes2$GOBP_Hsa_Inflammatory.response, "HP"),]
  rownames(NCC3)
  tNCC3 <- log(NCC3+1, 2)
  df <- tNCC3[!(rowMeans(tNCC3) < 2), ]
  dim(df)
  
  ijw <- calculate.correlation(df, doPerm = cor.perm,output.corTable = FALSE,output.permFDR = F)
  
  #### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
  run.par = doPar & (getDoParWorkers() == 1) 
  if (run.par) {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
  }
  
  ##### calculate PFN
  el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
  
  g <- graph.data.frame(el,directed = FALSE)
  
  ##### perform MCA clustering.
  MEGENA.output <- do.MEGENA(g,
                             mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = T,
                             min.size = 10,max.size = vcount(g)/2,
                             doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                             save.output = FALSE)
  
  ###### unregister cores as these are not needed anymore.
  if (getDoParWorkers() > 1) {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  
  summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                         mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                         min.size = 10,max.size = vcount(g)/2,
                                         annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                         output.sig = TRUE)
  
  if (!is.null(annot.table)) {
    # update annotation to map to gene symbols
    V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
    summary.output <- output[c("mapped.modules","module.table")]
    names(summary.output)[1] <- "modules"
  }
  
  print(head(summary.output$modules,3))
  summary.output$module.table$module.size
  
  #X11();
  pnet.obj <- plot_module(output.summary = summary.output,PFN = g,subset.module = "c1_4",
                          layout = "kamada.kawai",label.hubs.only = F,
                          gene.set = NULL,color.code =  "grey",
                          output.plot = FALSE,out.dir = "modulePlot", 
                          col.names = c("yellow", "red", "magenta", "green", "cyan"),label.scaleFactor = 20,
                          hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
  print(pnet.obj)
}  

