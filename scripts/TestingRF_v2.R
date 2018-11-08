source("~/Dropbox/codes/Harvard-BoeringerIngelheim/BI-compbio/R_codes/RFTyper/RF_functions.R")
#source("~/codes/BI-compbio/R_codes/RFTyper/RF_functions.R")
.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
library(tidyverse)
library(Seurat)

# PBMC + Pancreas Datasets testing for RF
load("~/Documents/RFTyper/pbmc3k_final.Rda")
load("~/Documents/RFTyper/pancreas.Robj")
#t-SNE plots showing cell types (with labels, etc) from pancreas and PBMC datasets
pdf("/Users/yasinkaymaz/Dropbox/HarvardInformatics/Tims-BI-Visit/Pancreas_Cells.tSNE.plot.pdf",width = 10,height = 8)
TSNEPlot(pancreas,group.by="assigned_cluster",do.label = T)
dev.off()

pdf("/Users/yasinkaymaz/Dropbox/HarvardInformatics/Tims-BI-Visit/PBMC_Cells.tSNE.plot.pdf",width = 10,height = 8)
TSNEPlot(pbmc,group.by="ClusterNames_0.6",do.label = T)
dev.off()

#model classification results (confusion matrix)
load("~/Documents/RFTyper/RF_model.pca.v3")
conf.mat <- rf$confusion %>% as.data.frame() %>% select(-class.error)
conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>%
  mutate(freq = 100*value/sum(value))
pdf("/Users/yasinkaymaz/Dropbox/HarvardInformatics/Tims-BI-Visit/Pancreas-PBMC_confusion-heatmap.pdf",width = 10,height = 8)
ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_y_discrete(name ="Predicted Cell Types")+
  scale_x_discrete(name ="Cell Types")
dev.off()


#Challenge of batch effects - uncorrected clustering of the 6 processed datasets showing batches
load("tasic.seurat.Robj")
load("romanov.seurat.Robj")
load("zeisel.seurat.Robj")
load("marques.seurat.Robj")
load("laManno.seurat.Robj")
load("gokce.seurat.Robj")


all6 <- MergeSeurat(object1=romanov, object2=zeisel, do.normalize =F)
all6 <- MergeSeurat(object1=all6, object2=marques, do.normalize =F)
all6 <- MergeSeurat(object1=all6, object2=laManno, do.normalize =F)
all6 <- MergeSeurat(object1=all6, object2=gokce, do.normalize =F)
all6 <- MergeSeurat(object1=all6, object2=tasic, do.normalize =F,do.scale = TRUE, do.center = TRUE)
all6 <- FindVariableGenes(all6, do.plot = F, display.progress = F)
hv.genes <- head(rownames(all6@hvg.info), 1000)
all6 <- RunPCA(all6, 
                        pc.genes = hv.genes, 
                        do.print = FALSE)
all6 <- FindClusters(all6, 
                              reduction.type = "pca", 
                              dims.use = 1:10, 
                              resolution = 1, 
                              print.output = FALSE, 
                              save.SNN = TRUE,
                              force.recalc = T)
all6 <- RunTSNE(all6, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)
Label="6_datasets_Before_Batch_Correction"
pdf(paste(Label,".plots.pdf",sep=""),width=8,height = 8)
PCAPlot(all6, dim.1 = 1, dim.2 = 2)
PCElbowPlot(all6)
TSNEPlot(all6, do.label = TRUE)
TSNEPlot(all6, do.label = TRUE,group.by ="tech")
dev.off()

# Experiments

setwd("/Users/yasinkaymaz/Documents/Harvard_Informatics/Data_Explore/mouse/FullSets")
load("mouseBrain.integrated.Aligned.seurat.Robj")

#First create subsets
study <- list("romanov", "laManno", "marques", "gokce", "tasic", "zeisel")
for (s in study){
  idlist <- mouseBrain.integrated@meta.data %>% as.tibble() %>% filter(tech==s) %>% select(orig.ident) %>% as.data.frame()
  print(head(idlist))
  #seu.obj <- paste(s,".obj",sep="");
  seu.obj <-  SubsetData(mouseBrain.integrated, cells.use = idlist$orig.ident, subset.raw=T,do.center = T,do.scale = T)
  print(head(seu.obj@meta.data))
  save(seu.obj, file=paste(s,".seurat.obj.subset.Robj",sep = ""))
  rm(seu.obj, idlist)
  gc()
}

rm(mouseBrain.integrated)


#Run from this: !!!!!!!!!!!!
source("~/codes/BI-compbio/R_codes/RFTyper/RF_functions.R")
.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
library(tidyverse)
library(Seurat)
library(randomForest)

args = commandArgs (T)
ExperimentName = args[1]


study <- list("romanov", "laManno", "marques", "gokce", "tasic", "zeisel")
for (s in study){
  load(paste(s,".seurat.obj.subset.Robj",sep = ""))
  assign(s, seu.obj)
  gc()
}

test.obj <- list(zeisel, tasic, laManno)
test.name <- c("zeisel", "tasic", "laManno")
newobjnames <- c("zeiselOut", "tasicOut", "laMannoOut")

#for (i in 1:length(newobjnames)){
i<-match(ExperimentName,newobjnames)
  print(i)
  print(test.obj[[i]])
  print(test.name[i])

# 1. Leave one out

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(
                list(tasic, romanov, laManno, marques, gokce),
                list(romanov, laManno, marques, gokce, zeisel),
                list(tasic, romanov, marques, gokce, zeisel)
                )

genes.use <- c()
for (x in 1:length(ob.list[[i]])) {
  ob.list[[i]][[x]] <- FindVariableGenes(ob.list[[i]][[x]], do.plot = F, display.progress = F)
  genes.use <- c(genes.use, head(rownames(ob.list[[i]][[x]]@hvg.info), 1000))
}

genes.use <- names(which(table(genes.use) > 1))
for (x in 1:length(ob.list[[i]])) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]][[x]]@scale.data)]
}


# Run multi-set CCA
ExperimentName <- RunMultiCCA(ob.list[[i]], genes.use = genes.use, num.ccs = 10)
save(ExperimentName, file=paste(test.name[i],"Out.integrated.seurat.Robj",sep=""))
# Run rare non-overlapping filtering
ExperimentName <- CalcVarExpRatio(object = ExperimentName, reduction.type = "pca",
                                         grouping.var = "tech", dims.use = 1:10)
ExperimentName <- SubsetData(ExperimentName, subset.name = "var.ratio.pca",
                                    accept.low = 0.5)

# Alignment
ExperimentName <- AlignSubspace(ExperimentName,
                                       reduction.type = "cca",
                                       grouping.var = "tech",
                                       dims.align = 1:10)

# t-SNE and Clustering
ExperimentName <- FindClusters(ExperimentName, 
                                      reduction.type = "cca.aligned",
                                      dims.use = 1:10, 
                                      save.SNN = T, 
                                      resolution = 1)

ExperimentName <- RunTSNE(ExperimentName,
                                 reduction.use = "cca.aligned",
                                 dims.use = 1:10,
                                 check_duplicates = FALSE)

save(ExperimentName, file=paste(test.name[i],"Out.integrated.seurat.Robj",sep=""))

# Visualization
pdf(paste(test.name[i],"Out.integrated.seurat.plots.pdf",sep=""))
MetageneBicorPlot(ExperimentName, grouping.var = "tech", dims.eval = 1:10)
TSNEPlot(ExperimentName, do.label = T)
TSNEPlot(ExperimentName, do.label = T, group.by ="tech")
dev.off()

# Model Training
sampleExpdata <- as.matrix(ExperimentName@data)
sampleClassdata <- ExperimentName@meta.data$res.1
rm(ExperimentName)

tmpname <- paste(test.name[i],"Out",sep = "")

prepareDataset(ExpressionData = sampleExpdata, CellLabels = sampleClassdata, do.splitTest = F, run.name = tmpname)
rm(sampleExpdata)
load(paste(tmpname,".trainingData.postPCA.data",sep=""))

CellTyperTrainer(trainingData.postPCA, tmpname)

load(paste(tmpname,"_rf.RF_model.Robj", sep = ""))




#on local
source("~/Dropbox/codes/Harvard-BoeringerIngelheim/BI-compbio/R_codes/RFTyper/RF_functions.R")
setwd("~/Documents/Harvard_Informatics/Data_Explore/mouse/FullSets/Experiments/")


#load("tasicOut.integrated.seurat.Robj")




load("tasicOut.RF_model.Robj")
rf <- modelname
rm(modelname)
load("tasic.seurat.obj.subset.Robj")
tasic <- seu.obj
rm(seu.obj)

CellLabels <- tasic@meta.data$res.1

testData <- as.data.frame(t(as.matrix(tasic@data)))
#It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
names(testData) <- make.names(names(testData))

indx <- sapply(testData, is.factor)
testData[indx] <- lapply(testData[indx], function(x) as.numeric(as.character(x)))
testData$CellType <- CellLabels
#Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
testData$CellType <- factor(testData$CellType)


testData <- droplevels(testData[,c(attributes(rf$terms)$term.labels, "CellType")])


CellTyper(testData, rf)
testSet <- testData
model <- rf

pred.test.prob <- predict(model, testSet, type = "prob")
pred.test.prob <- as.data.frame(pred.test.prob)

pred.test.out <- predict(model, testSet, type="response")

testPred <- pred.test.prob
testPred$Prior <- testSet$CellType

testPred$Prediction <- pred.test.out

confmat <- confusionMatrix(testPred$Prediction, testPred$Prior)
confmat$table

print(  
  confusionMatrix(data = testPred$Prediction,  
                  reference = testPred$Prior))

model$confusion[,"class.error"]

varImpPlot(model, sort = T, n.var=20,main="Top 10 - Variable Importance")








load("tasicOut.trainingData.postPCA.data")


load("tasic.seurat.obj.subset.Robj")

zeisel <- seu.obj
rm(seu.obj)
rf <- modelname
rm(modelname)

i=1
test.obj <- list(zeisel, tasic, laManno)
load("zeiselOut.trainingData.postPCA.data")
#colnames(trainingData.postPCA)
#xgene <- attributes(modelname$terms)$term.labels %>% gsub("X.", " ",.)

CellLabels <- zeisel@meta.data$res.1

testData <- as.data.frame(t(as.matrix(zeisel@data)))
#It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
names(testData) <- make.names(names(testData))

indx <- sapply(testData, is.factor)
testData[indx] <- lapply(testData[indx], function(x) as.numeric(as.character(x)))
testData$CellType <- CellLabels
#Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
testData$CellType <- factor(testData$CellType)


sort(colnames(trainingData.postPCA))

testData[1:10,colnames(trainingData.postPCA)]

testData[1:10,colnames(trainingData.postPCA)]
as.character(colnames(testData))
sapply(testData, is.factor)

TestExpdata <- 


rownames(test)

TestExpdata <- as.matrix(test.obj[[i]]@data)


prepareDataset(TestExpdata, TestClassdata, do.splitTest = T, percent.Testset=1.0, tmpname)

load(paste(tmpname,".testing.data",sep = ""))

modelname <- paste(tmpname, "_rf",sep = "")

testPred <- CellTyper(test, modelname)
testpred <- testPred %>% as.tibble() %>% group_by(Prior, Prediction) %>% summarize(n=n()) %>%
  mutate(freq = 100*n/sum(n))


data.class.err <- modelname$confusion %>% as.matrix() %>% as.data.frame() %>% select(class.error)
save(ExperimentName, file=paste(test.name[i],"Out.integrated.seurat.Robj",sep=""))

ExperimentName@meta.data$ClassPredError <-  data.class.err[ExperimentName@meta.data$res.1,]

errorSize <- as.data.frame(cbind(modelname$confusion[,"class.error"],
                                 table(trainingData.postPCA$CellType)))
colnames(errorSize) <- c("ClassError","ClassTrainingSize")
errorSize$CellTypeClass <- rownames(errorSize)


#Plots
pdf(paste(tmpname,".outputs.pdf",sep = ""),width=10,height = 8)

TSNEPlot(ExperimentName, do.label = TRUE)
TSNEPlot(ExperimentName, do.label = TRUE,group.by ="tech")


ExperimentName@meta.data %>% 
  as_tibble() %>% 
  ggplot(aes(as.integer(res.1), fill=tech)) + 
  geom_histogram(stat="count") + 
  scale_x_continuous(name ="Cluster IDs")+
  scale_y_continuous(name="Cell counts")

ExperimentName@meta.data %>% 
  as_tibble() %>% 
  ggplot(aes(tech, fill=res.1)) + 
  geom_histogram(stat="count") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(name="Cell counts")+
  scale_x_discrete(name ="Datasets")

ggplot(errorSize)  + 
  geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3")+
  geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
  scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize)))+
  theme(axis.text.x = element_text(angle = 0))+
  scale_x_discrete(name ="Cluster IDs")

ggplot(testpred, aes(Prior, Prediction, fill=freq)) +geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", 
                       name="% Predictions")

dev.off()


#}

