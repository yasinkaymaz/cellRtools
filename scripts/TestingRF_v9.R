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
}

zeisel.rank3.sub100 <- get(load("zeisel.rank3.sub100.seurat.Robj"))
# Training a Hierarchical RF model for 16 Rank-3 classes with equal number of cells, 886
zeisel.rank3.sub100 <- QuickSeurat(zeisel.rank3.sub100)#To add variable genes to var.genes
zeisel.rank3.sub100 <- SetAllIdent(zeisel.rank3.sub100, id = "TaxonomyRank3")#!!! important
zeisel.rank3.sub100 <- BuildClusterTree(zeisel.rank3.sub100, do.reorder = T, reorder.numeric = T)
CLoc.list <- AssessNodes(zeisel.rank3.sub100) #Creates hierarchical random forest classifiers
save(CLoc.list, file="CLoc.list.zeisel.rank3.sub100.Rdata")


#TESTING
zeisel.rank3.SM <- get(load("~/codes/test/RF/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))

zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))

zeisel5000.sm <- CellTyper(SeuratObject = zeisel5000subsampled, model = z.rank3.sm100, priorLabels = zeisel5000subsampled@meta.data$TaxonomyRank3, outputFilename = paste("zeisel5000subsampled",".predictions.w-SM100",sep = ""))
#Check accuracy...
#test with HRF 
zeisel5000.hrf <- HTyper2(SeuratObject = zeisel5000subsampled, models = CLoc.list)


