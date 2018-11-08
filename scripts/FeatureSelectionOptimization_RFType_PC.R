library(randomForest)
library(ggplot2)

#Functions:
selecteGenes.best.loadings <- function(p, pcs, num, caption="Highest loading"){
  #p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom
  load <- NULL
  topbotnames <- NULL
  for(i in 1:pcs){
    orderedpcai <- p$rotation[order(p$rotation[,i]),i]
    bot <- data.frame(genes=names(tail(orderedpcai,num)),bestgenes=tail(orderedpcai,num))
    top <- data.frame(genes=names(head(orderedpcai,num)),bestgenes=head(orderedpcai,num))
    load <- rbind(load, top, bot)
  }
  top <- tail(load[order(as.character(load$bestgenes)),],num/2)
  bot <- head(load[order(as.character(load$bestgenes)),],num/2)
  all <- rbind(top,bot)
  bestgenes <- droplevels(unique(all$genes))
  return(bestgenes)
}

#Second version:
selecteGenes.best.loadings <- function(p, pcs, num, caption="Highest loading"){
  #p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom
  #Weighted gene picking depending on PC number: Initial PCs give more genes.
  #For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  for(i in 1:pcs){
    orderedpcai <- p$rotation[order(abs(p$rotation[,i]),decreasing = TRUE),i]
    Gn <- round((pcs-i+1)*(num*2)/(pcs*(pcs+1)))
    TotalGenes = as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai,Gn)),bestgenes=head(orderedpcai,Gn))
    load <- rbind(load, top)
  }
  bestgenes <- unique(load$genes)
  return(bestgenes)
}

#PCrandomizer
PCrandomizer <- function(p, TotalPCs, pcs, num, caption="Highest loading"){
  #p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom
  #Weighted gene picking depending on PC number: Initial PCs give more genes. RANDOM PC - RANDOM GENES
  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  randomPCs <- sample(TotalPCs,pcs)
  for( i in randomPCs){
    orderedpcai <- p$rotation[,i]
    Gn <- round(i*num/sum(randomPCs))
    TotalGenes = as.numeric(TotalGenes) + Gn
    sampledGenes <- sample(orderedpcai,Gn)
    top <- data.frame(genes=names(sampledGenes),bestgenes=sampledGenes)
    load <- rbind(load, top)
  }
  bestgenes <- unique(load$genes)
  
  return(bestgenes)
}


#MODEL TRAINING
#************************************#
#load training dataset
setwd("~/codes/test")


load("trainingData.data")
#Filter three low accuracy cell types out from data
trainingData <- trainingData[!(trainingData$CellType %in% c("t_cell", "schwann", "epsilon")),]

tsub <- trainingData[,!(names(trainingData) %in% c("CellType"))]
tsub <- droplevels(tsub)
indx <- sapply(tsub, is.factor)
tsub[indx] <- lapply(tsub[indx], function(x) as.numeric(as.character(x)))
head(tsub[1:10,1:10])

#First split the data into training and test sets
smp_size <- floor(0.80 * nrow(tsub))
## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(tsub)), size = smp_size)
train <- tsub[train_ind, ]
test <- tsub[-train_ind, ]
# Then do PCA on only training set:
train <- train[ , apply(train, 2, var) != 0]



#Optimize PCs and number of genes to use as features
load("PCA.train.data")

errorsTable <- NULL
AccuracyTable <- NULL
RuntimeTable <- NULL
rand.errorsTable <- NULL
rand.AccuracyTable <- NULL

for (classnum in seq(from=1, to=20, by=1)){
	FeatNum=1000
	#CLASSNumber=19
	print(classnum)
	start_time <- Sys.time()

	pca.genes <- as.character(selecteGenes.best.loadings(pcatrain,classnum,FeatNum))

	trainfiltered <- train[, pca.genes]

	trainingData.postPCA <- data.frame(trainfiltered, CellType=trainingData[rownames(trainfiltered),"CellType"])
	trainingData.postPCA <- droplevels(trainingData.postPCA)

	rf = randomForest(CellType~., data = trainingData.postPCA, norm.votes = TRUE,importance=TRUE, proximity = TRUE, ntree=500)
	save(rf, file=paste("RF_model.pca.test_",classnum,"_PCs.PCA.Robj",sep = ""))

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

	acc.table <- data.frame(PCsIncluded=classnum, Accuracy=acc)
	AccuracyTable <- rbind(AccuracyTable, acc.table)

	errorSize <- data.frame(PCsIncluded=classnum,ClassError = rf$confusion[,"class.error"])
	errorSize$CellTypeClass <- rownames(errorSize)

	errorsTable <- rbind(errorsTable, errorSize)

	end_time <- Sys.time()
	runtime <- end_time - start_time
	runt <- data.frame(PCsIncluded=classnum,Runtime = runtime)

	RuntimeTable <- rbind(RuntimeTable, runt)
	
	#Random Gene Selection
	rand.genes <- as.character(PCrandomizer(pcatrain,19,classnum,FeatNum ))
	trainfiltered <- train[, rand.genes]
	trainingData.postPCA <- data.frame(trainfiltered, CellType=trainingData[rownames(trainfiltered),"CellType"])
	trainingData.postPCA <- droplevels(trainingData.postPCA)
	rf = randomForest(CellType~., data = trainingData.postPCA, norm.votes = TRUE,importance=TRUE, proximity = TRUE, ntree=500)
	testSet <- test[,rand.genes]
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
  
	rand.acc <- 100*length(SuccessPredictions[,1]) / (length(SuccessPredictions[,1]) +length(FailedPredictions[,1]) )
  	rand.acc.table <- data.frame(rand.FeatureNumber=length(rand.genes), rand.Accuracy=rand.acc)
  	rand.AccuracyTable <- rbind(rand.AccuracyTable, rand.acc.table)
  
  	rand.errorSize <- data.frame(rand.FeatureNumber=length(rand.genes),rand.ClassError = rf$confusion[,"class.error"])
  	rand.errorSize$rand.CellTypeClass <- rownames(rand.errorSize)
  	rand.errorsTable <- rbind(rand.errorsTable, rand.errorSize)
  
  	save(AccuracyTable, file="AccuracyTable_PCsTest.Robj")
  	save(errorsTable, file="ErrorsTable_PCsTest.Robj")
  	save(RuntimeTable, file="RuntimeTable_PCsTest.Robj")
  
  	save(rand.AccuracyTable, file="rand.AccuracyTable_PCsTest.Robj")
	save(rand.errorsTable, file="rand.ErrorsTable_PCsTest.Robj")

}

AccuracyTable <- cbind(AccuracyTable, rand.AccuracyTable)
errorsTable <- cbind(errorsTable, rand.errorsTable)

save(AccuracyTable, file="AccuracyTable_PCsTest.Robj")
save(errorsTable, file="ErrorsTable_PCsTest.Robj")
save(RuntimeTable, file="RuntimeTable_PCsTest.Robj")

pacc <- ggplot(AccuracyTable)+geom_line(aes(x=PCsIncluded, y=Accuracy),color="blue")+ylim(c(0,100))+geom_line(aes(x=PCsIncluded, y=rand.Accuracy),color="red")
perr <- ggplot(errorsTable, aes(x=PCsIncluded, y=ClassError, color=CellTypeClass))+geom_line()+ylim(c(0,1))+geom_point(aes(shape=CellTypeClass))

ggsave(plot=pacc,width = 12,height = 8, dpi=200, filename="PCs_accuracyPlot.pdf", useDingbats=FALSE )
ggsave(plot=perr, width = 20,height = 8, dpi=200, filename="PCs_ClassErrorPlot.pdf", useDingbats=FALSE )


