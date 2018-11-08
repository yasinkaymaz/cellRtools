
source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")
library(tidyverse)
library(Seurat)
setwd("~/Documents/RFTyper/test3/")
#Load multiple PBMC datasets

dir='~/Documents/RFTyper/Datasets/'

datas <- c(
  "pbmc4K_CR.2.1.0",
  "pbmc8K_CR.2.1.0"
)

for (d in datas){
  print(d)
  exp.data <- Read10X(data.dir = paste(dir,d,"/filtered_gene_bc_matrices/GRCh38/",sep = ""))
  SeuratWrapper(d, exp.data, Normalize = T, Label = d)
  rm(exp.data)
}

load("../pbmc3k_final.Rda")



# Model Training
sampleExpdata <- as.matrix(pbmc@data)
sampleClassdata <- pbmc@meta.data$ClusterNames_0.6

run.name = "pbmc3K"
prepareDataset(ExpressionData = sampleExpdata, CellLabels = sampleClassdata, do.splitTest = F, run.name = run.name)

rm(sampleExpdata)
load(paste(run.name,".trainingData.postPCA.data",sep=""))

CellTyperTrainer(trainingData.postPCA, run.name)

load(paste(run.name,".RF_model.Robj", sep = ""))
rf <- modelname
rm(modelname)
head(pbmc@meta.data)

library(rfUtilities)
dim(trainingData.postPCA)


rf$confusion
TSNEPlot(pbmc,group.by="res.0.6")
FeaturePlot(pbmc,features.plot = "nGene")
FeaturePlot(pbmc,features.plot = "nUMI")
FeaturePlot(pbmc,features.plot = "percent.mito")


TestExpData <- t(as.matrix(pbmc8K_CR.2.1.0@data))

head(TestExpData[1:10,1:10])
CellTyper(TestExpData, rf)
#example.lab <- as.data.frame(PredictionOutput[,"Prediction",drop=F])

source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")


TestExpData <- t(as.matrix(pbmc8K_CR.2.1.0@data))

CellTyper(TestExpData, rf)

pbmc8K_CR.2.1.0@meta.data$Prediction <- PredictionOutput$Prediction

pdf("pbmc8k.predictions.tsne.pdf",width =10,height = 12)
TSNEPlot(pbmc, group.by="ClusterNames_0.6",do.label=T)
PCAPlot(pbmc, group.by="ClusterNames_0.6",do.label=T)
TSNEPlot(pbmc8K_CR.2.1.0, group.by="Prediction",do.label=T)

FeaturePlot(object = pbmc8K_CR.2.1.0, 
            features.plot = c("MS4A1", "GNLY", "IL7R", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
varImpPlot(rf, sort = T, n.var=50,main="Top 10 - Variable Importance")

dev.off()


#Try picking DE genes for features 
pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

top100.markers <- pbmc.markers %>% group_by(cluster) %>% top_n(100, avg_logFC)

# Model Training
sampleExpdata <- as.matrix(pbmc@data)
sampleClassdata <- pbmc@meta.data$ClusterNames_0.6

run.name = "pbmc3Kdge"
training.genes <- unique(top100.markers$gene)

prepareDataset(ExpressionData = sampleExpdata, CellLabels = sampleClassdata, do.splitTest = F, run.name = run.name, featureGeneSet = training.genes)

rm(sampleExpdata)
load(paste(run.name,".trainingData.postPCA.data",sep=""))

CellTyperTrainer(trainingData.postPCA, run.name)

load(paste(run.name,".RF_model.Robj", sep = ""))

rf.dge <- modelname

CellTyper(TestExpData, rf.dge)

pbmc8K_CR.2.1.0@meta.data$Prediction <- PredictionOutput$Prediction

pdf("pbmc8kdge.predictions.tsne.pdf",width =10,height = 12)
TSNEPlot(pbmc, group.by="ClusterNames_0.6",do.label=T)
PCAPlot(pbmc, group.by="ClusterNames_0.6",do.label=T)
TSNEPlot(pbmc8K_CR.2.1.0, group.by="Prediction",do.label=T)

FeaturePlot(object = pbmc8K_CR.2.1.0, 
            features.plot = c("MS4A1", "GNLY", "IL7R", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
varImpPlot(rf.dge, sort = T, n.var=50,main="Top 10 - Variable Importance")

dev.off()


pbmc8K_CR.2.1.0 <- RunTSNE(pbmc8K_CR.2.1.0, dims.use = 1:10, do.fast = TRUE,perplexity=40, check_duplicates = FALSE)

pdf("pbmc8k-tsne.predictions.tsne.pdf",width =10,height = 12)

TSNEPlot(pbmc8K_CR.2.1.0, group.by="Prediction",do.label=T)

dev.off()



source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")

bestScoreExtractor(PredictionOutput)

pbmc8K_CR.2.1.0@meta.data$BestScores <- bestscoresList$BestVotesPercent

pdf("pbmc8k-tsne.predictions.tsne.pdf",width =10,height = 12)
TSNEPlot(pbmc8K_CR.2.1.0, group.by="Prediction",do.label=T)
FeaturePlot(pbmc8K_CR.2.1.0, features.plot = "BestScores")

ggplot(data=pbmc8K_CR.2.1.0@meta.data,aes(y=BestScores,x=Prediction))+
  geom_violin()

dev.off()




#plot the relationship between the class sizes in the training set and the distribution of best scores per class from the predictions
table(trainingData.postPCA[,"CellType"]) %>% as.data.frame()

library(reshape2)
PredictionOutput %>% as.tibble() %>% filter(.,Prediction == "CD4 T cells") %>% select(-Prediction) %>% filter(`CD4 T cells` <= 0.5)


pdf("Diagnosis.predictions.1.pdf",width =10,height = 12)
PredictionOutput %>% as.tibble() %>% filter(.,Prediction == "CD4 T cells") %>% select(-Prediction) %>% filter(`CD4 T cells` <= 0.5) %>% melt() %>% ggplot(., aes(y=value, x=variable)) + geom_violin()
dev.off()


output <- cbind(pbmc8K_CR.2.1.0@meta.data[,c("nGene", "nUMI")], PredictionOutput)
head(output)

pdf("Diagnosis.predictions.2.pdf",width =10,height = 12)
output %>% as.data.frame() %>% filter(Prediction == "CD4 T cells" & `CD4 T cells` <= 0.5) %>% ggplot(., aes(x=nGene, y=`B cells`))+geom_point()
dev.off()


output <- cbind(pbmc8K_CR.2.1.0@meta.data[,c("nGene", "nUMI", "BestScores")], PredictionOutput)
head(output)

pdf("Diagnosis.predictions.3.pdf",width =10,height = 12)
output %>% as.data.frame() %>% ggplot(., aes(x=nUMI, y=BestScores))+geom_point()
dev.off()



important.genes <- rf$importance %>% 
  as.tibble() %>% 
  add_column(Genes=rownames(rf$importance)) %>% 
  arrange(., desc(MeanDecreaseGini)) %>% 
  top_n(n = 200,wt=MeanDecreaseGini) %>% 
  select(Genes) %>% 
  as.list()

important.genes <- gsub("\\.", "-", important.genes$Genes)
important.genes <- important.genes[c(which(! important.genes %in% c("RP11-290F20-3","TMEM66")))]



#Modify Seurat object meta.data slot with predictions
source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")

testpredictionoutput <- CellTyper(TestExpData, rf.dge)

#pbmc8K_CR.2.1.0@meta.data$Prediction <- PredictionOutput$Prediction

bestscoresList <- bestScoreExtractor(testpredictionoutput)

pbmc8K_CR.2.1.0@meta.data$BestScores <- bestscoresList$BestVotesPercent

output <- cbind(pbmc8K_CR.2.1.0@meta.data[,c("nGene", "nUMI", "BestScores")], testpredictionoutput)
head(output)
pbmc8K_CR.2.1.0@meta.data <- output
pdf("Diagnosis.predictions.5.pdf",width= 16,height = 16)
FeaturePlot(object = pbmc8K_CR.2.1.0, 
            features.plot = rf$classes, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
TSNEPlot(pbmc8K_CR.2.1.0, group.by="Prediction",do.label=T)
FeaturePlot(pbmc8K_CR.2.1.0, features.plot = "BestScores")

dev.off()


#Testing 4k pbmc data

pbmc4K_CR.2.1.0 <- CellTyper(SeuratObject = pbmc4K_CR.2.1.0, model = rf)

PlotPredictions(SeuratObject = pbmc4K_CR.2.1.0, model = rf, save.pdf = T,outputFilename = "PBMC4k.Predictions")

#Test with Pancreas dataset
load("~/Documents/RFTyper/Datasets/pancreas.Robj")
head(pancreas@meta.data)
pancreas <- CellTyper(SeuratObject = pancreas, model=rf)                       
PlotPredictions(SeuratObject = pancreas, model = rf, save.pdf = T,outputFilename = "pancreas.Predictions")

#### STRaTA           ####
#Class size balancing

source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")
run.name="rf.with.strata"
CellTyperTrainer(trainingData.postPCA, run.name)

load(paste(run.name,".RF_model.Robj", sep = ""))

rf.strata <- modelname

head(pbmc4K_CR.2.1.0@meta.data)

pbmc4K_CR.2.1.0 <- CellTyper(SeuratObject = pbmc4K_CR.2.1.0, model = rf.strata)

PlotPredictions(SeuratObject = pbmc4K_CR.2.1.0, model = rf, save.pdf = T,outputFilename = paste(run.name,"PBMC4k.Predictions",sep = ""))

#merge two datasets
ob.list <- list(pbmc, pbmc4K_CR.2.1.0, pbmc8K_CR.2.1.0)
ob.list <- list(pbmc, pbmc4K_CR.2.1.0)
pbmc7K.integrated <- SeuratCCAmerger(listofObjects = ob.list)

# Visualization
TSNEPlot(pbmc7K.integrated, do.label = T)
TSNEPlot(pbmc7K.integrated, do.label = T, group.by ="dataSource")

head(pbmc7K.integrated@meta.data)
pbmc7K.integrated <- CellTyper(SeuratObject = pbmc7K.integrated, model = rf.strata)

PlotPredictions(SeuratObject = pbmc7K.integrated, model = rf.strata, outputFilename = "PBMC7k.Predictions")

datas <- c(
  "pbmc33K_CR.1.1.0"
)

for (d in datas){
  print(d)
  exp.data <- Read10X(data.dir = paste(dir,d,"/filtered_gene_bc_matrices/hg19/",sep = ""))
  SeuratWrapper(d, exp.data, Normalize = T, Label = d)
  rm(exp.data)
}

pbmc6K_CR.1.1.0@meta.data
TSNEPlot(pbmc6K_CR.1.1.0, do.label = T)

pbmc6K_CR.1.1.0 <- CellTyper(SeuratObject = pbmc6K_CR.1.1.0, model = rf.strata)

PlotPredictions(SeuratObject = pbmc6K_CR.1.1.0, model = rf.strata, outputFilename = "pbmc6K_CR.1.1.0.Predictions")














library(randomForest)
library(AUC)

make.data = function(N=1000) {
  X = data.frame(replicate(6,rnorm(N))) #six features
  y = X[,1]^2+sin(X[,2]) + rnorm(N)*1 #some hidden data structure to learn
  rare.class.prevalence = 0.1
  y.class = factor(y<quantile(y,c(rare.class.prevalence))) #10% TRUE, 90% FALSE
  return(data.frame(X,y=y.class))
}

#make some data structure
train.data = make.data()

#1 - Balancing by voting rule, AUC of ROC will be unchanged...
rare.class.prevalence = 0.1
rf.cutoff = randomForest(y~.,data=train.data,cutoff=c(1-rare.class.prevalence,rare.class.prevalence))
print(rf.cutoff)

#2 - Balancing by sampling stratification
nRareSamples = 1000 * rare.class.prevalence
rf.strata = randomForest(y~.,data=train.data,strata=train.data$y,
                         sampsize=c(nRareSamples,nRareSamples))
print(rf.strata)

#3 - Balancing by class-weight during training.
rf.classwt = randomForest(y~.,data=train.data,classwt=c(0.0005,1000))
print(rf.classwt)

#view OOB-CV specificity and sensitiviy
plot(roc(rf.cutoff$votes[,2],train.data$y),main="black default, red stata, green classwt")
plot(roc(rf.strata$votes[,2],train.data$y),col=2,add=T)
plot(roc(rf.classwt$votes[,2],train.data$y),col=3,add=T)


#make test.data and remove random sample until both classes are equally prevalent
test.data = make.data(N=50000)
test.data.balanced = test.data[-sample(which(test.data$y=="FALSE"))[1:40000],]

#print prediction performance %predicted correct:
sapply(c("rf.cutoff","rf.strata","rf.classwt"),function(a.model) {
  mean(test.data.balanced$y == predict(get(a.model), newdata=test.data.balanced))
})






