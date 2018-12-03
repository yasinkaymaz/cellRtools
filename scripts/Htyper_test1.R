.libPaths("~/biotools/Rlibs")
source("~/codes/cellRtools/main/functions.R")
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyverse)

zeisel.rank1.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank1.rf.RF_model.Robj"))
zeisel.rank2.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank2.rf.RF_model.Robj"))
zeisel.rank3.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank3.rf.RF_model.Robj"))
zeisel.rank4.rf <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Zeisel2018/zeisel.rank4.rf.RF_model.Robj"))

#zeisel.rank1.rfcv <- get(load("~/codes/test/RF/zeisel.rank1.rfcv.RF_model_notImproved.Robj"))
#zeisel.rank2.rfcv <- get(load("~/codes/test/RF/zeisel.rank2.rfcv.RF_model_notImproved.Robj"))
#zeisel.rank3.rfcv <- get(load("~/codes/test/RF/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
#zeisel.rank4.rf <- get(load("~/codes/test/RF/redo_rank4/zeisel.rank4.rfcv.RF_model_notImproved.Robj"))

models.list <- list(zeisel.rank1.rf,zeisel.rank2.rf,zeisel.rank3.rf,zeisel.rank4.rf)

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
