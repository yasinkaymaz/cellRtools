
source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")
library(tidyverse)
library(Seurat)
setwd("~/Documents/RFTyper/test3/")
#Load multiple PBMC datasets

dir='~/Documents/RFTyper/Datasets/'

datas <- c(
"CD14posMonocytes",
"CD19posBCells",
"CD34posCells",
"CD4posHelperTCells",
"CD4posCD25posRegulatoryTCells",
"CD4posCD45RAposCD25negNaiveTcells",
"CD4posCD45ROposMemoryTCells",
"CD56posNaturalKillerCells",
"CD8posCytotoxicTcells",
"CD8posCD45RAposNaiveCytotoxicTCells"
)

for (d in datas){
  print(d)
  exp.data <- Read10X(data.dir = paste(dir,d,"/",sep = ""))
  SeuratWrapper(d, exp.data, Normalize = T, Label = d)
  rm(exp.data)
}


merged <- ScaleData(merged,genes.use = hv.genes)

merged <- RunPCA(merged,
                 pc.genes = hv.genes,
                 do.print = FALSE)

merged <- FindClusters(merged,
                       reduction.type = "pca",
                       dims.use = 1:10,
                       resolution = 1,
                       print.output = FALSE,
                       save.SNN = TRUE,
                       force.recalc = T)
merged <- RunTSNE(merged, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)

head(merged@meta.data)

load("rf.with.strata.RF_model.Robj")
rf <- modelname

merged <- CellTyper(SeuratObject = merged, model = rf)
PlotPredictions(SeuratObject = merged, model = rf, outputFilename = "merged.CD14-CD19-CD34.predictions")

TSNEPlot(merged, do.label = T)


require(scales)

# Create vector with levels of object@ident
merged <- SetAllIdent(object = merged, id = "Prediction")
identities <- levels(merged@ident)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# Plot the tSNE plot with the default ggplot2 colors
TSNEPlot(object = merged, do.return = T) + 
  scale_color_manual(values = my_color_palette)

FeatureHeatmap(object = merged, features.plot = "BestVotesPercent", group.by = "Prediction", sep.scale = F, pt.size = 1.5, pch.use = 20)


#optimize rf model
load("pbmc3K.trainingData.postPCA.data")

rf.improved <- CellTyperTrainer(trainingData = trainingData.postPCA, run.name = "improved.RF")


merged <- CellTyper(SeuratObject = merged, model = rf.improved)
PlotPredictions(SeuratObject = merged, model = rf.improved, outputFilename = "rf.improved.merged.CD14-CD19-CD34.predictions")

merged@meta.data$KS.UniformTest.p <- apply(merged@meta.data[,5:12], 1, function(x) ks.test(as.numeric(x),"punif",0,1)$p.value )
FeatureHeatmap(object = merged, features.plot = "KS.UniformTest.p", group.by = "Prediction", sep.scale = F, pt.size = 1.5, pch.use = 20)
library(entropy)
merged@meta.data$entropy <- apply(merged@meta.data[,5:12], 1, function(x) entropy(as.numeric(x)) )
plot(merged@meta.data$BestVotesPercent, merged@meta.data$entropy,pch=20)

head(merged@meta.data)
FeatureHeatmap(object = merged, features.plot = "entropy", group.by = "res.1", sep.scale = F, pt.size = 1.5, pch.use = 20)


merged@meta.data$KLe <- apply(merged@meta.data[,5:12], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/8,8)) )
plot(merged@meta.data$BestVotesPercent, merged@meta.data$KLe,pch=20)

merged@meta.data[which(merged@meta.data$KLe < 0.5),]
hist(merged@meta.data$KLe)

merged@meta.data %>% as.tibble() %>% mutate(Prediction2 = ifelse( KLe <= 0.5, "Undetermined", "Detected")) %>% as.data.frame()
merged@meta.data %>% as.tibble() %>% select(Prediction) %>% c() %>% table()
merged@meta.data %>% as.tibble() %>% filter(KLe < 0.5) %>% select(Prediction) %>% c() %>% table()

#Run with updated function
merged <- CellTyper(SeuratObject = merged, model = rf.improved)
head(merged@meta.data)
FeatureHeatmap(object = merged, features.plot = "KLe", group.by = "PredictionStatus", sep.scale = F, pt.size = 1.5, pch.use = 20)
FeaturePlot(merged,features.plot = "KLe",no.legend = F,cols.use = c("gray","red"))

FeaturePlot(merged,features.plot = "BestVotesPercent",no.legend = F,cols.use = c("gray","red"))

FeatureHeatmap(object = merged, features.plot = "KLe", group.by = "PredictionStatus", sep.scale = F, pt.size = 1.5, pch.use = 20, cols.use = c("gray","purple"))
FeaturePlot(merged,features.plot = "KLe",no.legend = F,cols.use = c("gray","purple"))

PlotPredictions(merged,model = rf.improved,outputFilename = "rf.improved.merged.CD14-CD19-CD34.predictions.wKle")



load("improved.RF.RF_model.Robj")
rf.improved <- rf

datas <- c(
  "CD14posMonocytes",
  "CD19posBCells",
  "CD34posCells",
  "CD4posHelperTCells",
  "CD4posCD25posRegulatoryTCells",
  "CD4posCD45RAposCD25negNaiveTcells",
  "CD4posCD45ROposMemoryTCells",
  "CD56posNaturalKillerCells",
  "CD8posCytotoxicTcells",
  "CD8posCD45RAposNaiveCytotoxicTCells"
)

for (d in datas){
  print(d)
  outfile=paste("rf.improved",d,"predictions",sep = ".")
  exp.data <- Read10X(data.dir = paste(dir,d,"/",sep = ""))
  d <- SeuratWrapper(SeuratObjName = d, ExpData = exp.data, Normalize = T, Label = d)
  rm(exp.data)
  d <- CellTyper(SeuratObject = d, model = rf.improved)
  PlotPredictions(SeuratObject = d, model = rf.improved, outputFilename = outfile)
  rm(d)
}

load("pbmc3K.trainingData.postPCA.data")
rf.improved


source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")

rf.imp.cved <- CellTyperTrainer(trainingData = trainingData.postPCA,run.name = "rf.imp.cved",improve = F )

cvresults <- rfcv(trainx = trainingData.postPCA[,-c(length(trainingData.postPCA[1,]))], trainingData.postPCA$CellType, cv.fold = 5)

rf <- randomForest(x = trainingData.postPCA[,-c(length(trainingData.postPCA[1,]))], y = trainingData.postPCA$CellType, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData.postPCA$CellType)))

rf.cv <- rf.crossValidation(x = rf, xdata = trainingData.postPCA[,-c(length(trainingData.postPCA[1,]))],norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData.postPCA$CellType)), bootstrap = T)

#https://panglaodb.se/samples.html?species=mouse&protocol=all%20protocols&sort=mostrecent
load("~/Downloads/SRA617020_SRS2566665.sparse.RData")
rownames(sm) <- make.unique(str_split(rownames(sm),"_",simplify = T)[,1])

sm <- as.matrix(sm)
dim(sm)

bcell.sm <- SeuratWrapper(SeuratObjName = "sm", ExpData = sm, Label = "sm",Normalize = T)
dim(bcell.sm@meta.data)
bcell.sm <- CellTyper(SeuratObject = bcell.sm, model = rf.improved)
PlotPredictions(SeuratObject = bcell.sm,model = rf.improved,outputFilename = "Bcells.fromPanglo.predictions")


load("~/Documents/RFTyper/Datasets/SRA.final/SRA701877_SRS3279684.sparse.RData")
rownames(sm) <- make.unique(str_split(rownames(sm),"_",simplify = T)[,1])

sm <- as.matrix(sm)
dim(sm)
panccells <- SeuratWrapper(SeuratObjName = "sm", ExpData = sm, Label = "panc", Normalize = T)
head(panccells@meta.data)
panccells <- CellTyper(SeuratObject = panccells, model = rf.improved)
PlotPredictions(SeuratObject = panccells,model = rf.improved,outputFilename = "PancreaticCells-1.fromPanglo.predictions")

panccells@meta.data %>% as.tibble() %>% filter(PredictionStatus == "Detected") %>% as.data.frame()
FeatureHeatmap(object = panccells, features.plot = "KLe", group.by = "PredictionStatus", sep.scale = F, pt.size = 1.5, pch.use = 20, cols.use = c("gray","purple"))

plot(panccells@meta.data$BestVotesPercent, panccells@meta.data$KLe,pch=20)
rf.improved

datas <- c("pbmc33K_CR.1.1.0")
for (d in datas){
  print(d)
  outfile=paste("rf.improved",d,"predictions",sep = ".")
  exp.data <- Read10X(data.dir = paste(dir,d,"/",sep = ""))
  d <- SeuratWrapper(SeuratObjName = d, ExpData = exp.data, Normalize = T, Label = d)
  rm(exp.data)
  d <- CellTyper(SeuratObject = d, model = rf.improved)
  PlotPredictions(SeuratObject = d, model = rf.improved, outputFilename = outfile)
  rm(d)
}

d@meta.data


d.15 <- CellTyper(SeuratObject = d.15, model = rf.improved)
PlotPredictions(SeuratObject = d.15, model = rf.improved, outputFilename = "d.15.predictions")

d.15@meta.data %>% as.tibble() %>% filter(PredictionStatus == "Detected") %>% as.data.frame()

d.15@meta.data %>% as.tibble() %>% select(-c(nGene, nUMI, orig.ident, res.1, KLe, Prediction, PredictionStatus, BestVotesPercent)) %>% apply(., 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])


d.x <- SubsetData(d,subset.name = "res.1",accept.value="8",subset.raw = T)

source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")
d.x <- CellTyper(SeuratObject = d.x, model = rf.improved)
plot(d.x@meta.data$Diff, d.x@meta.data$KLe, pch=20)
d.x@meta.data %>% as.tibble() %>% filter(PredictionStatus == "Undetermined") %>% as.data.frame()
plot(d.x@meta.data$Diff, d.x@meta.data$nUMI, pch=20)
plot(d.x@meta.data$Diff, d.x@meta.data$nGene, pch=20)
hist(d.x@meta.data$KLe)
hist(d.x@meta.data$Diff)
PlotPredictions(SeuratObject = d.x, model = rf.improved, outputFilename = "d.8.predictions")
d.x@meta.data %>% as.tibble() %>% filter(Prediction != "NK cells") %>% as.data.frame()

FeatureHeatmap(object = d.x, features.plot = "Diff", group.by = "PredictionStatus", sep.scale = F, pt.size = 1.5, pch.use = 20, cols.use = c("gray","purple"))

FeaturePlot(d.15,features.plot = "Diff",no.legend = F,cols.use = c("gray","red"))

class_n=8
d.x@meta.data %>% as.tibble() %>% select(-c(nGene, nUMI, orig.ident, res.1)) %>% mutate(., Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class_n), "Unclassified", Prediction))  %>% as.data.frame()

