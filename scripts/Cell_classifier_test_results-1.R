library(randomForest)
#library(Seurat)
#library(caret)
library(ggplot2)

#Functions:



#MODEL TRAINING
#************************************#
#load training dataset
setwd("~/Documents/RFTyper/ModelOptimization/version2/")

intervals <- as.numeric(c("10","20","30","40","50","60","70","80","90","100","110","120","130","140","150","160","170","180","190","200","300","500","700","900","1100","1300","1500","1700","1900","2100","2300","2500","2700","2900","3100","3300","3500","3700","3900","4100","4300","4500","4700","4900"))
errorsTable <- NULL
AccuracyTable <- NULL
RuntimeTable <- NULL

for (FeatNum in intervals){
  print(FeatNum)
  start_time <- Sys.time()
  
  CLASSNumber=19
  #pca.genes <- as.character(selecteGenes.best.loadings(pcatrain,CLASSNumber,FeatNum))
  load(paste("RF_model.pca.test_",FeatNum,"_Genes.PCA.Robj",sep = ""))
  print("Loaded")
  pca.genes <- rownames(rf$importance)
  
  testSet <- test[,pca.genes]
  testSet <- data.frame(testSet, CellType=trainingData[rownames(testSet),"CellType"])
  pred.test.prob <- predict(rf, testSet, type = "prob")
  pred.test.prob <- as.data.frame(pred.test.prob)
  pred.test.out <- predict(rf, test, type="response")
  testPred <- pred.test.prob
  testPred$Prior <- testSet$CellType
  testPred$Prediction <- pred.test.out
  testPred$Prior <- droplevels(testPred$Prior)
  FailedPredictions <- testPred[as.character(testPred$Prediction) != as.character(testPred$Prior),]
  SuccessPredictions <- testPred[as.character(testPred$Prediction) == as.character(testPred$Prior),]
  
  acc <- 100*length(SuccessPredictions[,1]) / (length(SuccessPredictions[,1]) +length(FailedPredictions[,1]) )
  
  acc.table <- data.frame(FeatureNumber=length(pca.genes), Accuracy=acc)
  AccuracyTable <- rbind(AccuracyTable, acc.table)
  
  errorSize <- data.frame(FeatureNumber=length(pca.genes),ClassError = rf$confusion[,"class.error"])
  errorSize$CellTypeClass <- rownames(errorSize)
  
  errorsTable <- rbind(errorsTable, errorSize)
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  runt <- data.frame(FeatureNumber=length(pca.genes),Runtime = runtime)
  
  RuntimeTable <- rbind(RuntimeTable, runt)
}

save(AccuracyTable, file="AccuracyTable_FeatureTest.Robj")
save(errorsTable, file="ErrorsTable_FeatureTest.Robj")
save(RuntimeTable, file="RuntimeTable_FeatureTest.Robj")

load("~/Documents/RFTyper/ModelOptimization/version2/AccuracyTable_FeatureTest.Robj")
load("~/Documents/RFTyper/ModelOptimization/version2/ErrorsTable_FeatureTest.Robj")
load("~/Documents/RFTyper/ModelOptimization/version2/RuntimeTable_FeatureTest.Robj")

load("~/Documents/RFTyper/ModelOptimization/version2/rand.AccuracyTable_FeatureTest.Robj")
load("~/Documents/RFTyper/ModelOptimization/version2/rand.ErrorsTable_FeatureTest.Robj")

AccuracyTable <- cbind(AccuracyTable, rand.AccuracyTable)
ggplot(AccuracyTable)+geom_line(aes(x=FeatureNumber, y=Accuracy),color="blue")+ylim(c(0,100))+geom_line(aes(x=FeatureNumber, y=rand.Accuracy),color="red")

load("~/Documents/RFTyper/ModelOptimization/version2/AccuracyTable_PCsTest.Robj")
load("~/Documents/RFTyper/ModelOptimization/version2/rand.AccuracyTable_PCsTest.Robj")

AccuracyTable <- cbind(AccuracyTable, rand.AccuracyTable)
ggplot(AccuracyTable)+geom_line(aes(x=PCsIncluded, y=Accuracy),color="blue")+ylim(c(0,100))+geom_line(aes(x=PCsIncluded, y=rand.Accuracy),color="red")



RuntimeTable[10:25,]$Runtime <- 60*RuntimeTable[10:25,]$Runtime





pacc <- ggplot(AccuracyTable, aes(x=FeatureNumber, y=Accuracy))+geom_line(stat="identity")+ylim(c(75,100))
perr <- ggplot(errorsTable, aes(x=FeatureNumber, y=ClassError, color=CellTypeClass))+geom_line()+ylim(c(0,1))+geom_point(aes(shape=CellTypeClass))
prunt <- ggplot(RuntimeTable, aes(x=FeatureNumber, y=Runtime))+geom_line(stat="identity")
ggplot(RuntimeTable, aes(x=FeatureNumber, y=Runtime))+geom_point(stat="identity")

ggsave(plot=pacc,width = 12,height = 8, dpi=200, filename="accuracyPlot.pdf", useDingbats=FALSE )
ggsave(plot=perr, width = 20,height = 8, dpi=200, filename="ClassErrorPlot.pdf", useDingbats=FALSE )
ggsave(plot=prunt, width = 20,height = 8, dpi=200, filename="PredictionRunTime_vs_FeatureNumPlot.pdf", useDingbats=FALSE )



#PC number optimization:
load("~/Documents/RFTyper/ModelOptimization/version2/AccuracyTable_PCsTest.Robj")
load("~/Documents/RFTyper/ModelOptimization/version2/ErrorsTable_PCsTest.Robj")
load("~/Documents/RFTyper/ModelOptimization/version2/RuntimeTable_PCsTest.Robj")

pacc <- ggplot(AccuracyTable, aes(x=PCsIncluded, y=Accuracy))+geom_line(stat="identity")+ylim(c(75,100))
perr <- ggplot(errorsTable, aes(x=PCsIncluded, y=ClassError, color=CellTypeClass))+geom_line()+ylim(c(0,1))+geom_point(aes(shape=CellTypeClass))
prunt <- ggplot(RuntimeTable, aes(x=PCsIncluded, y=Runtime))+geom_line(stat="identity")



ggsave(plot=pacc,width = 12,height = 8, dpi=200, filename="PCs_accuracyPlot.pdf", useDingbats=FALSE )
ggsave(plot=perr, width = 20,height = 8, dpi=200, filename="PCs_ClassErrorPlot.pdf", useDingbats=FALSE )


tdata <- read.delim("/Users/yasinkaymaz/Documents/whatsup/sorted.text.count.txt",header = 1)
hist(tdata$Number_of_Texts_to_GrooupChat,breaks = 500,col = "red",lwd=2)
