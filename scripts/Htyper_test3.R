.libPaths("~/biotools/Rlibs")
source("~/codes/cellRtools/main/functions.R")
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyverse)
library(DropSeq.util)


zeisel.rank1.rfcv <- get(load("~/codes/test/RF/zeisel.rank1.rfcv.RF_model_notImproved.Robj"))
zeisel.rank2.rfcv <- get(load("~/codes/test/RF/zeisel.rank2.rfcv.RF_model_notImproved.Robj"))
zeisel.rank3.rfcv <- get(load("~/codes/test/RF/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
zeisel.rank4.rfcv <- get(load("~/codes/test/RF/redo_rank4/rank4.cv/zeisel.rank4.rfcv.RF_model_notImproved.Robj"))

TaxonomyRankTable <- read.delim("~/codes/test/RF/TaxonomyRank.tree.zeisel.txt", header=F)

models.list <- list(zeisel.rank1.rfcv,zeisel.rank2.rfcv,zeisel.rank3.rfcv,zeisel.rank4.rfcv)

ho <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Ho2018/Ho2018.seurat.Robj"))
ho <- HTyper22(SeuratObject = ho, taxTable=TaxonomyRankTable, models = models.list, priorLabels = ho@meta.data$res.1, outputFilename = paste("ho",".rfcv.predictions.w-HM",sep = ""))
rm(ho)

tasic2018 <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Tasic2018/Tasic2018.seurat.Robj"))
tasic2018 <- HTyper22(SeuratObject = tasic2018, taxTable=TaxonomyRankTable, models = models.list, priorLabels = tasic2018@meta.data$res.1, outputFilename = paste("tasic2018",".rfcv.predictions.w-HM",sep = ""))
rm(tasic2018)

Scer <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Saunders/P60Cerebellum_ALT.seurat.Robj"))
Scer <- HTyper22(SeuratObject = Scer, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Scer@meta.data$res.1, outputFilename = paste("Saunders-P60Cerebellum_ALT",".rfcv.predictions.w-HM",sep = ""))
rm(Scer)

Scor <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Saunders/P60Cortex_noRep5_FRONTALonly.seurat.Robj"))
Scor <- HTyper22(SeuratObject = Scor, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Scor@meta.data$res.1, outputFilename = paste("Saunders-P60Cortex_noRep5_FRONTALonly",".rfcv.predictions.w-HM",sep = ""))
rm(Scor)

Sstr <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Saunders/P60Striatum.seurat.Robj"))
Sstr <- HTyper22(SeuratObject = Sstr, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Sstr@meta.data$res.1, outputFilename = paste("Saunders-P60Striatum",".rfcv.predictions.w-HM",sep = ""))
rm(Sstr)

Ssn <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Saunders/P60SubstantiaNigra.seurat.Robj"))
Ssn <- HTyper22(SeuratObject = Ssn, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Ssn@meta.data$res.1, outputFilename = paste("Saunders-P60SubstantiaNigra",".rfcv.predictions.w-HM",sep = ""))
rm(Ssn)

Sth <- get(load("/n/home13/yasinkaymaz/LabSpace/testdata/GEO/Saunders/P60Thalamus.seurat.Robj"))
Sth <- HTyper22(SeuratObject = Sth, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Sth@meta.data$res.1, outputFilename = paste("Saunders-P60Thalamus",".rfcv.predictions.w-HM",sep = ""))
rm(Sth)
