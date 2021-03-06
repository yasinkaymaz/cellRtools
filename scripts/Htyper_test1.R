.libPaths("~/biotools/Rlibs")
source("~/codes/cellRtools/main/functions.R")
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyverse)
library(DropSeq.util)

zeisel.rank1.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank1.rf.RF_model.Robj"))
zeisel.rank2.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank2.rf.RF_model.Robj"))
zeisel.rank3.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank3.rf.RF_model.Robj"))
zeisel.rank4.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank4.rf.RF_model.Robj"))

#zeisel.rank1.rfcv <- get(load("~/codes/test/RF/zeisel.rank1.rfcv.RF_model_notImproved.Robj"))
#zeisel.rank2.rfcv <- get(load("~/codes/test/RF/zeisel.rank2.rfcv.RF_model_notImproved.Robj"))
#zeisel.rank3.rfcv <- get(load("~/codes/test/RF/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
#zeisel.rank4.rf <- get(load("~/codes/test/RF/redo_rank4/rank4.cv/zeisel.rank4.rfcv.RF_model_notImproved.Robj"))

models.list <- list(zeisel.rank1.rf,zeisel.rank2.rf,zeisel.rank3.rf,zeisel.rank4.rf)



zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))
zeisel5000subsampled <- HTyper(SeuratObject = zeisel5000subsampled, models = models.list, priorLabels = zeisel5000subsampled@meta.data$TaxonomyRank4, outputFilename = paste("zeisel5000subsampled",".rf.predictions.w-HM",sep = ""))
rm(zeisel5000subsampled)

Astrocytes <- get(load("~/LabSpace/testdata/GEO/Saunders/Astrocytes.seurat.Robj"))
Astrocytes <- HTyper(SeuratObject = Astrocytes, models = models.list, priorLabels = Astrocytes@meta.data$Region, outputFilename = paste("Astrocytes",".rf.predictions.w-HM",sep = ""))
rm(Astrocytes)

Endothelial <- get(load("~/LabSpace/testdata/GEO/Saunders/Endothelial.seurat.Robj"))
Endothelial <- HTyper(SeuratObject = Endothelial, models = models.list, priorLabels = Endothelial@meta.data$Region, outputFilename = paste("Endothelial",".rf.predictions.w-HM",sep = ""))
rm(Endothelial)

FibroblastLike <- get(load("~/LabSpace/testdata/GEO/Saunders/FibroblastLike.seurat.Robj"))
FibroblastLike <- HTyper(SeuratObject = FibroblastLike, models = models.list, priorLabels = FibroblastLike@meta.data$Region, outputFilename = paste("FibroblastLike",".rf.predictions.w-HM",sep = ""))
rm(FibroblastLike)

Microglia <- get(load("~/LabSpace/testdata/GEO/Saunders/Microglia.seurat.Robj"))
Microglia <- HTyper(SeuratObject = Microglia, models = models.list, priorLabels = Microglia@meta.data$Region, outputFilename = paste("Microglia",".rf.predictions.w-HM",sep = ""))
rm(Microglia)

Mural <- get(load("~/LabSpace/testdata/GEO/Saunders/Mural.seurat.Robj"))
Mural <- HTyper(SeuratObject = Mural, models = models.list, priorLabels = Mural@meta.data$Region, outputFilename = paste("Mural",".rf.predictions.w-HM",sep = ""))
rm(Mural)

Oligodendrocytes <- get(load("~/LabSpace/testdata/GEO/Saunders/Oligodendrocytes.seurat.Robj"))
Oligodendrocytes <- HTyper(SeuratObject = Oligodendrocytes, models = models.list, priorLabels = Oligodendrocytes@meta.data$Region, outputFilename = paste("Oligodendrocytes",".rf.predictions.w-HM",sep = ""))
rm(Oligodendrocytes)

Polydendrocytes <- get(load("~/LabSpace/testdata/GEO/Saunders/Polydendrocytes.seurat.Robj"))
Polydendrocytes <- HTyper(SeuratObject = Polydendrocytes, models = models.list, priorLabels = Polydendrocytes@meta.data$Region, outputFilename = paste("Polydendrocytes",".rf.predictions.w-HM",sep = ""))
rm(Polydendrocytes)

# Testing models on real mixture datasets
#Testing on negative control
MCA_Muscle <- get(load("~/LabSpace/testdata/GEO/Han2018/MCA_Muscle.seurat.Robj"))
MCA_Muscle <- HTyper(SeuratObject = MCA_Muscle, models = models.list, priorLabels = MCA_Muscle@meta.data$Annotation, outputFilename = paste("MCA_Muscle",".rf.predictions.w-HM",sep = ""))
rm(MCA_Muscle)

Gokce2016 <- get(load("~/LabSpace/testdata/GEO/Gokce2016/Gokce2016.seurat.Robj"))
Gokce2016 <- HTyper(SeuratObject = Gokce2016, models = models.list, priorLabels = Gokce2016@meta.data$type, outputFilename = "Gokce2016.rf.predictions.w-HM")
rm(Gokce2016)

Harris2018 <- get(load("~/LabSpace/testdata/GEO/Harris2018/Harris2018_corrected.seurat.Robj"))
Harris2018 <- HTyper(SeuratObject=Harris2018, priorLabels=Harris2018@meta.data$Comments, models=models.list, outputFilename="Harris2018.rf.predictions.w-HM")
rm(Harris2018)

Romanov2017 <- get(load("~/LabSpace/testdata/GEO/Romanov2017/Romanov2017.seurat.Robj"))
Romanov2017 <- HTyper(SeuratObject=Romanov2017, models = models.list, priorLabels = Romanov2017@meta.data$level1.class, outputFilename="Romanov2017.rf.predictions.w-HM")
rm(Romanov2017)

Hrvatin <- get(load("~/LabSpace/testdata/GEO/Hrvatin2018/Hrvatin.seurat.Robj"))
Hrvatin <- HTyper(SeuratObject=Hrvatin, models = models.list, priorLabels = Hrvatin@meta.data$celltype, outputFilename="Hrvatin.rf.predictions.w-HM")
rm(Hrvatin)

Hochgerner2018 <- get(load("~/LabSpace/testdata/GEO/Hochgerner2018/Hochgerner2018.seurat.Robj"))
Hochgerner2018 <- HTyper(SeuratObject=Hochgerner2018, models = models.list, priorLabels = Hochgerner2018@meta.data$res.1, outputFilename="Hochgerner2018.rf.predictions.w-HM")
rm(Hochgerner2018)

Campbell2017  <- get(load("~/LabSpace/testdata/GEO/Campbell2017/Campbell2017.seurat.Robj"))
Campbell2017 <- HTyper(SeuratObject = Campbell2017, models = models.list, priorLabels = Campbell2017@meta.data$res.1, outputFilename = "Cambell2017.rf.predictions.w-HM")
rm(Campbell2017)

Tasic2016  <- get(load("~/LabSpace/testdata/GEO/Tasic2016/Tasic2016.seurat.Robj"))
Tasic2016 <- HTyper(SeuratObject=Tasic2016, models = models.list, priorLabels = Tasic2016@meta.data$res.1,  outputFilename="Tasic2016.rf.predictions.w-HM")
rm(Tasic2016)


hdatafiles <- list.files(pattern = "*htable.Rdata")
library(tidyverse)
 library(alluvial)
 library(ggalluvial)

pdf("DataSets.prediction-crosscheck.pdf",width = 14, height = 10)
for(i in 1:length(hdatafiles)){
  type <- strsplit(hdatafiles[i],split = "\\.")[[1]][1]
  print(type)
  load(hdatafiles[i])
  p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2], axis3 = names(crx)[3], axis4 = names(crx)[4], axis5 = names(crx)[5], axis5 = names(crx)[6] )) +
      geom_alluvium(aes_string(fill = names(crx)[6]), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = names(crx), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")
  print(p5)
}
dev.off()



zmode1 <- zeisel.rank1.rfcv

#Recursive training/Prediction
zeisel.rank1.sub <- get(load("zeisel.rank1.sub.seurat.Robj"))

imp <- as.data.frame(zeisel.rank1.rfcv$finalModel$importance)
features <- rownames(head(imp[order(imp$MeanDecreaseGini, decreasing=T),],200))
cells <- rownames(zeisel.rank1.rfcv$finalModel$votes)

R1.LR2 <- zeisel.rank1.sub@meta.data[cells,c("TaxonomyRank1","TaxonomyRank2")]
colnames(R1.LR2) <- c("R1","CellType")

ExpressionData=as.data.frame(as.matrix(zeisel.rank1.sub@data))

trainingData <- as.data.frame(t(ExpressionData))
  #It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
  names(trainingData) <- make.names(names(trainingData))
  trainingData <- trainingData[,features]
  trainingData <- cbind(trainingData, R1.LR2)
  #Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
  trainingData$CellType <- factor(trainingData$CellType)

zmode2 <- RecursiveTrainer(trainingData = trainingData, run.name="z.modeL2")
