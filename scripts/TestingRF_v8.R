.libPaths("~/biotools/Rlibs")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/test-functions.R")

setwd("/n/home13/yasinkaymaz/codes/test/RF/HRF_test8")

#The Goal is to compare single model prediction accuracy to hierarchical-randomForest model predictions.

zeisel.rank3.sub <- get(load("~/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank3.sub.seurat.Robj"))
run='pass'

# Training a single RF model for 16 Rank-3 classes with equal number of cells, 886
if(run == 'pass'){
  print("not running this...")
}else{
  #TODO: add code how processed Zeisel data upto this point.
  #Load Zeisel data:
  load("Zeisel2018.seurat.Robj")#check directory
  #Subset cells to balance Rank3 classes:
  cells <- NULL; classes <- names(table(zeisel@meta.data$TaxonomyRank3)); min_n <- min(table(zeisel@meta.data$TaxonomyRank3))
  for(i in 1:length(classes)){
    cells <- c(cells,sample(rownames(zeisel@meta.data[which(zeisel@meta.data$TaxonomyRank3 == classes[i]),]),min_n))
  }
  zeisel.rank3.sub <- SubsetData(object = zeisel, cells.use = cells,do.clean=T)
  #Train a single RF model with 16 Classes:
  z.rank3.rfcv <- CellTyperTrainer2(ExpressionData = zeisel.rank3.sub@data, CellLabels = zeisel.rank3.sub@meta.data$TaxonomyRank3, run.name = "zeisel.rank3.rfcv", do.splitTest = F, PCs = 40, improve.rf = T)
}

zeisel.rank3.rfcv <- get(load("~/codes/test/RF/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))


#TESTING
zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))
zeisel5000subsampled <- HTyper22(SeuratObject = zeisel5000subsampled, taxTable=TaxonomyRankTable, models = models.list, priorLabels = zeisel5000subsampled@meta.data$TaxonomyRank4, outputFilename = paste("zeisel5000subsampled",".rfcv.predictions.w-HM",sep = ""))
rm(zeisel5000subsampled)




