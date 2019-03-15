.libPaths("~/biotools/Rlibs")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/test-functions.R")

setwd("/n/home13/yasinkaymaz/codes/test/RF/HRF_test9")

#The Goal is to compare single model prediction accuracy to hierarchical-randomForest model predictions.

run='pass'

# Training a single RF model for 16 Rank-3 classes with equal number of cells, 886
if(run == 'pass'){
  print("not running this...")
}else{
  #TODO: add code how processed Zeisel data upto this point.
  #Load the zeisel subset data already present:
  zeisel.rank3.sub <- get(load("~/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank3.sub.seurat.Robj"))
  #Subset cells to balance Rank3 classes:
  cells <- NULL; classes <- names(table(zeisel.rank3.sub@meta.data$TaxonomyRank3)); min_n <- 100;
  for(i in 1:length(classes)){
    cells <- c(cells,sample(rownames(zeisel.rank3.sub@meta.data[which(zeisel.rank3.sub@meta.data$TaxonomyRank3 == classes[i]),]),min_n))
  }
  zeisel.rank3.sub100 <- SubsetData(object = zeisel.rank3.sub, cells.use = cells,do.clean=T)
  save(zeisel.rank3.sub100, file="zeisel.rank3.sub100.seurat.Robj")
  #Train a single RF model with 16 Classes:
  z.rank3.sm100 <- CellTyperTrainer2(ExpressionData = zeisel.rank3.sub100@data, CellLabels = zeisel.rank3.sub100@meta.data$TaxonomyRank3, run.name = "zeisel.rank3.sm100", do.splitTest = F, PCs = 40, improve.rf = F)
  save(z.rank3.sm100, file="z.rank3.sm100.Rmod")
  # Training a Hierarchical RF model for 16 Rank-3 classes with equal number of cells, 886
  zeisel.rank3.sub100 <- QuickSeurat(zeisel.rank3.sub100)#To add variable genes to var.genes
  zeisel.rank3.sub100 <- SetAllIdent(zeisel.rank3.sub100, id = "TaxonomyRank3")#!!! important
  #zeisel.rank3.sub100 <- BuildClusterTree(zeisel.rank3.sub100, do.reorder = T, reorder.numeric = T)
  zeisel.rank3.sub100 <- BuildClusterTree(zeisel.rank3.sub100, do.reorder = T)
  CLoc.list <- AssessNodes(zeisel.rank3.sub100) #Creates hierarchical random forest classifiers
  save(CLoc.list, file="CLoc.list.zeisel.rank3.sub100.Rdata")
  
}

zeisel.rank3.sub100 <- get(load("zeisel.rank3.sub100.seurat.Robj"))

ztraintest <- BuildClusterTree(zeisel.rank3.sub100, do.reorder = T)


#TESTING
zeisel.rank3.SM <- get(load("z.rank3.sm100.Rmod"))

zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))

zeisel5000.sm <- CellTyper2(SeuratObject = zeisel5000subsampled, model = z.rank3.sm100, priorLabels = zeisel5000subsampled@meta.data$TaxonomyRank3, outputFilename = paste("zeisel5000subsampled",".predictions.w-SM100",sep = ""))
#Check accuracy...
table(zeisel5000.sm@meta.data[,c("Prior", "Intermediate")])
#test with HRF 
zeisel5000.hrf <- HTyper2(SeuratObject = zeisel5000subsampled, models = CLoc.list, tree=zeisel.rank3.sub100@cluster.tree[[1]])
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

pdf("conf.mat.pdf",width = 12,height = 8)
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

confmat <- table(zeisel5000.hrf@meta.data[,c("TaxonomyRank3", "HRFPrediction")]  %>% mutate(HRFPrediction = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(HRFPrediction)] ))
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

confmat <- table(zeisel5000.hrf@meta.data[,c("TaxonomyRank3", "Decision")]  %>% mutate(Decision = ztraintest@cluster.tree[[1]]$tip.label[as.numeric(Decision)] ))
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
dev.off()

pdf("tree.pdf", width = 15, height = 15)
zeisel.rank3.sub100 <- SetAllIdent(zeisel.rank3.sub100, id = "TaxonomyRank3")#!!! important
zeisel.rank3.sub100 <- BuildClusterTree(zeisel.rank3.sub100, do.reorder = T)
dev.off()

