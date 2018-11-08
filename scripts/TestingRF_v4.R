.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
source("/n/home13/yasinkaymaz/codes/scTyper/RFTyper/RF_functions.R")
library(tidyverse)
library(Seurat)
setwd("/n/home13/yasinkaymaz/codes/test/brainCellMarkers/")
#Load multiple PBMC datasets

dir='/n/home13/yasinkaymaz/codes/test/brainCellMarkers/Datasets/'

datasets <- c(
  "pbmc4K_CR.2.1.0",
  "pbmc8K_CR.2.1.0",
  "pbmc33K_CR.1.1.0",
  "pbmc3K_CR.1.1.0",
  "pbmc68K_CR.1.0.0",
  "pbmc6K_CR.1.1.0",
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

for (d in datasets){
  print(d)
  exp.data <- Read10X(data.dir = paste(dir,d,"/",sep = ""))
  SeuratWrapper(d, exp.data, Normalize = T, Label = d)
  rm(exp.data)
}



# Model Training
load(paste(dir,"/pbmc3k_final.Rda",sep = ""))
sampleExpdata <- as.matrix(pbmc@data)
sampleClassdata <- pbmc@meta.data$ClusterNames_0.6

run.name = "pbmc3K"
prepareDataset(ExpressionData = sampleExpdata, CellLabels = sampleClassdata, do.splitTest = F, run.name = run.name)

rm(sampleExpdata)
load(paste(run.name,".trainingData.postPCA.data",sep=""))

CellTyperTrainer(trainingData.postPCA, run.name)

load(paste(run.name,".RF_model.Robj", sep = ""))
rf.strata <- modelname
rm(modelname)
head(pbmc@meta.data)


# Testing
ob.list <- list(pbmc4K_CR.2.1.0, pbmc8K_CR.2.1.0, pbmc33K_CR.1.1.0, pbmc3K_CR.1.1.0, pbmc68K_CR.1.0.0, pbmc6K_CR.1.1.0)

pbmc122K.integrated <- SeuratCCAmerger(listofObjects = ob.list)

pbmc122K.integrated <- CellTyper(SeuratObject = pbmc122K.integrated, model = rf.strata)

PlotPredictions(SeuratObject = pbmc122K.integrated, model = rf.strata, outputFilename = "PBMC122k.Predictions")

save(pbmc122K.integrated, file="pbmc122K.integrated.seurat.Robj")

rm(ob.list, sampleClassdata, trainingData.postPCA, )

gc()

oblist <- list(pbmc4K_CR.2.1.0, pbmc8K_CR.2.1.0, pbmc33K_CR.1.1.0, pbmc3K_CR.1.1.0, pbmc68K_CR.1.0.0, pbmc6K_CR.1.1.0, CD14posMonocytes, CD19posBCells, CD34posCells, CD4posHelperTCells, CD4posCD25posRegulatoryTCells, CD4posCD45RAposCD25negNaiveTcells, 
               CD4posCD45ROposMemoryTCells, CD56posNaturalKillerCells, CD8posCytotoxicTcells, CD8posCD45RAposNaiveCytotoxicTCells)

for(i in 1:length(oblist)){
  if(i == 1){
    merged <- oblist[[i]]
  }else{
    merged <- MergeSeurat(object1 = merged, object2 = oblist[[i]], add.cell.id1 = "merged1", add.cell.id2 = "merged2")
  }
}

rm(oblist, pbmc4K_CR.2.1.0, pbmc8K_CR.2.1.0, pbmc33K_CR.1.1.0, pbmc3K_CR.1.1.0, pbmc68K_CR.1.0.0, pbmc6K_CR.1.1.0, CD14posMonocytes, CD19posBCells, CD34posCells, CD4posHelperTCells, CD4posCD25posRegulatoryTCells, CD4posCD45RAposCD25negNaiveTcells,
               CD4posCD45ROposMemoryTCells, CD56posNaturalKillerCells, CD8posCytotoxicTcells, CD8posCD45RAposNaiveCytotoxicTCells )

gc()

merged <- FindVariableGenes(merged, do.plot = F, display.progress = F)
hv.genes <- head(rownames(merged@hvg.info), 1000)
merged <- ScaleData(merged)

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
save(merged, file="merged.seurat.objects.Robj")


