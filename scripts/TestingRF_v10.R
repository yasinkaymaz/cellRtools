.libPaths("~/biotools/Rlibs")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/test-functions.R")

setwd("/n/home13/yasinkaymaz/codes/test/RF/HRF_test10")

#The Goal is to compare single model prediction accuracy to hierarchical-randomForest model predictions.

run='pass'

# Training a single RF model for 39 Rank-3 classes with equal number of cells, 47
if(run == 'pass'){
  print("not running this...")
}else{
  #TODO: add code how processed Zeisel data upto this point.
  #Load the zeisel subset data already present:
  zeisel.rank4.sub <- get(load("/n/home13/yasinkaymaz/codes/test/RF/redo_rank4/zeisel.rank4.sub.Rdata"))
  #Subset cells to balance Rank3 classes: 
  cells <- NULL; classes <- names(table(zeisel.rank4.sub@meta.data$TaxonomyRank4)); min_n <- min(table(zeisel.rank4.sub@meta.data$TaxonomyRank4));
  for(i in 1:length(classes)){
    cells <- c(cells,sample(rownames(zeisel.rank4.sub@meta.data[which(zeisel.rank4.sub@meta.data$TaxonomyRank4 == classes[i]),]),min_n))
  }
  zeisel.rank4.sub47 <- SubsetData(object = zeisel.rank4.sub, cells.use = cells,do.clean=T)
  save(zeisel.rank4.sub47, file="zeisel.rank4.sub47.seurat.Robj")
  #Train a single RF model with 16 Classes:
  z.rank4.sm47 <- CellTyperTrainer2(ExpressionData = zeisel.rank4.sub47@data, CellLabels = zeisel.rank4.sub47@meta.data$TaxonomyRank4, run.name = "zeisel.rank4.sm47", do.splitTest = F, PCs = 40, improve.rf = F)
  save(z.rank4.sm47, file="z.rank4.sm47.Rmod")
  
  # Training a Hierarchical RF model for 16 Rank-3 classes with equal number of cells, 886
  zeisel.rank4.sub47 <- QuickSeurat(zeisel.rank4.sub47)#To add variable genes to var.genes
  zeisel.rank4.sub47 <- SetAllIdent(zeisel.rank4.sub47, id = "TaxonomyRank4")#!!! important
  zeisel.rank4.sub47 <- BuildClusterTree(zeisel.rank4.sub47, do.reorder = T, reorder.numeric = T)
  #zeisel.rank4.sub47 <- BuildClusterTree(zeisel.rank4.sub47, do.reorder = T)
  CLoc.list <- AssessNodes(zeisel.rank4.sub47) #Creates hierarchical random forest classifiers
  save(CLoc.list, file="CLoc.list.zeisel.rank4.sub47.Rdata")
  
}


pdf("tree.pdf", width = 15, height = 20)
zeisel.rank4.sub47 <- SetAllIdent(zeisel.rank4.sub47, id = "TaxonomyRank4")
ztraintest <- BuildClusterTree(zeisel.rank4.sub47, do.reorder = T)
dev.off()

#TESTING
z.rank4.sm47 <- get(load("z.rank4.sm47.Rmod"))

zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))

zeisel5000.sm <- CellTyper2(SeuratObject = zeisel5000subsampled, model = z.rank4.sm47, priorLabels = zeisel5000subsampled@meta.data$TaxonomyRank4, outputFilename = paste("zeisel5000subsampled",".predictions.w-SM47",sep = ""))
#Check accuracy...
table(zeisel5000.sm@meta.data[,c("Prior", "Intermediate")])
#test with HRF
zeisel5000.hrf <- HTyper2(SeuratObject = zeisel5000subsampled, models = CLoc.list, tree=zeisel.rank4.sub47@cluster.tree[[1]])
table(zeisel5000.hrf@meta.data[,c("HRFPrediction", "Decision")])
table(zeisel5000.hrf@meta.data[,c("TaxonomyRank3", "Decision")])

zeisel5000.hrf@meta.data[,c("HRFPrediction", "Decision", "TaxonomyRank3")] %>%  mutate(Decision = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(Decision)] ) %>% mutate(HRFPrediction = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(HRFPrediction)] )



tipBias <- function(tree, confmat){
  #confmat is table() of Prior/Prediction comparison
  err.rate <- NULL
  for(i in 1:length(tree$tip.label)){
    nn <- length(GetAncestorsPath(tree=tree, i)[[2]])
    leafErr <- 1-confmat[tree$tip.label[i],tree$tip.label[i]]/sum(confmat[tree$tip.label[i],])
    print(paste(i,nn,leafErr,sep = "    "))
    err.rate <- c(err.rate, leafErr)
  }
  err.rate <- data.frame(err.rate)
  rownames(err.rate) <- tree$tip.label
  return(err.rate)
}
#tipBias(tree=ztraintest@cluster.tree[[1]], confmat=confmat)

zeisel.rank4.sub47 <- SetAllIdent(zeisel.rank4.sub47, id = "TaxonomyRank4")
ztraintest <- BuildClusterTree(zeisel.rank4.sub47, do.reorder = T)


pdf("conf.mat.pdf",width = 16,height = 12)
confmat <- table(zeisel5000.sm@meta.data[,c("Prior", "Intermediate")])
err <- tipBias(tree=ztraintest@cluster.tree[[1]], confmat=confmat)
pheatmap::pheatmap(confmat,
                   cluster_cols = F,
                   cluster_rows = F,
                   main = paste("Single Model Prediction. Cohen's Kappa =",round(DescTools::CohenKappa(confmat),digits = 3),sep = " "),
                   scale = "none",
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f",
                   annotation_col = err)

in.sm.pred <- rownames(zeisel5000.sm@meta.data[which(zeisel5000.sm@meta.data$Prior != zeisel5000.sm@meta.data$Intermediate),])
in.mat <- table(zeisel5000.sm@meta.data[which(rownames(zeisel5000.sm@meta.data) %in% in.sm.pred),c("Prior", "Intermediate")]  )
pheatmap::pheatmap(in.mat,
                   cluster_cols = F,
                   cluster_rows = F,
                   main = paste("Single Model Prediction. Incorrect predictions =",length(in.sm.pred),sep = " "),
                   scale = "none",
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f")

confmat <- table(zeisel5000.hrf@meta.data[,c("TaxonomyRank4", "HRFPrediction")]  %>% mutate(HRFPrediction = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(HRFPrediction)] ))
err <- tipBias(tree=ztraintest@cluster.tree[[1]], confmat=confmat)
pheatmap::pheatmap(confmat,
                   cluster_cols = F,
                   cluster_rows = F,
                   main = paste("HRF Prediction. Cohen's Kappa =",round(DescTools::CohenKappa(confmat),digits = 3),sep = " "),
                   scale = "none",
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f",
                   annotation_col = err)

in.hrf.pred <- zeisel5000.hrf@meta.data[,c("cells","TaxonomyRank4", "HRFPrediction")] %>%
  mutate(HRFPrediction = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(HRFPrediction)] ) %>%
  filter(TaxonomyRank4 != HRFPrediction ) %>%
  select(cells)
in.mat <- table(zeisel5000.hrf@meta.data[which(rownames(zeisel5000.hrf@meta.data) %in% in.hrf.pred$cells),c("TaxonomyRank4", "HRFPrediction")]  %>% mutate(HRFPrediction = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(HRFPrediction)] ))
pheatmap::pheatmap(in.mat,
                   cluster_cols = F,
                   cluster_rows = F,
                   main = paste("HRF Prediction. Incorrect predictions =",length(in.hrf.pred$cells),sep = " "),
                   scale = "none",
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f")

confmat <- table(zeisel5000.hrf@meta.data[,c("TaxonomyRank4", "Decision")]  %>% mutate(Decision = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(Decision)] ))
err <- tipBias(tree=ztraintest@cluster.tree[[1]], confmat=confmat)
pheatmap::pheatmap(confmat, 
                   cluster_cols = F,
                   cluster_rows = F,
                   main = paste("Gated Decision. Cohen's Kappa =",round(DescTools::CohenKappa(confmat),digits = 3),sep = " "),
                   scale = "none",
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f",
                   annotation_col = err)
in.hrf.pred <- zeisel5000.hrf@meta.data[,c("cells","TaxonomyRank4", "Decision")] %>%
  mutate(Decision = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(Decision)] ) %>%
  filter(TaxonomyRank4 != Decision ) %>%
  select(cells)
in.mat <- table(zeisel5000.hrf@meta.data[which(rownames(zeisel5000.hrf@meta.data) %in% in.hrf.pred$cells),c("TaxonomyRank4", "Decision")]  %>% mutate(Decision = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(Decision)] ))
pheatmap::pheatmap(in.mat,
                   cluster_cols = F,
                   cluster_rows = F,
                   main = paste("Gated Decision. Incorrect predictions =",length(in.hrf.pred$cells),sep = " "),
                   scale = "none",
                   display_numbers = T,
                   number_color = "black",
                   number_format = "%.0f")

dev.off()


in.hrf.pred <- zeisel5000.hrf@meta.data[,c("cells","TaxonomyRank4", "HRFPrediction")] %>%
  mutate(HRFPrediction = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(HRFPrediction)] ) %>%
  filter(TaxonomyRank4 != HRFPrediction ) %>%
  select(cells)
length(in.hrf.pred$cells)


