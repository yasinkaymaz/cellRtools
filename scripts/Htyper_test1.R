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
#zeisel.rank4.rf <- get(load("~/codes/test/RF/redo_rank4/zeisel.rank4.rfcv.RF_model_notImproved.Robj"))

models.list <- list(zeisel.rank1.rf,zeisel.rank2.rf,zeisel.rank3.rf,zeisel.rank4.rf)


setwd("~/LabSpace/testdata/GEO/Saunders/")
#Saunders pure populations for testing Positive controls
# datafiles <- list.files(pattern="H_1stRound_CrossTissue_.*txt.gz$")
# metafiles <- list.files(pattern="H_1stRound_CrossTissue_.*.RDS$")
#
# for(i in 1:length(datafiles)){
#   region <- strsplit(datafiles[i],split = "_")[[1]][4]
#   print(region)
#   data <- loadSparseDge(datafiles[i])
#   print(dim(data))
#   meta <- readRDS(metafiles[i])
#   print(dim(meta))
#   meta$Region <- region
#   print(head(meta))
#   seuratobj <- SeuratWrapper(ExpData = data, ProjectLabel = region, NewMeta = meta, Normalize = T, scale.only.var = T, PCs = 20, dump.files = F)
#   seuratobj <- HTyper(SeuratObject = seuratobj, models = models.list, priorLabels = seuratobj@meta.data$Region, outputFilename = paste(region,".rf.predictions.w-HM",sep = ""))
#   remove(seuratobj)
# }

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

pdf("PositiveControls.prediction-crosscheck.pdf",width = 14, height = 10)
for(i in 1:length(hdatafiles)){
  type <- strsplit(hdatafiles[i],split = "\\.")[[1]][1]
  print(type)
  load(hdatafiles[i])
  crx.f <- crx %>% mutate(freq = n*100 / sum(n)) %>% filter(freq > 1)
  p5 <- ggplot(crx.f,aes_string(y = "n", axis1 = names(crx.f)[1], axis2 = names(crx.f)[2], axis3 = names(crx.f)[3], axis4 = names(crx.f)[4], axis5 = names(crx.f)[5] )) +
    geom_alluvium(aes_string(fill = names(crx.f)[5]), width = 0, knot.pos = 1/4) +
    guides(fill = FALSE)+
    geom_stratum(width = 1/12, fill = "grey", color = "red") +
    geom_label(stat = "stratum", label.strata = TRUE ) +
    ylab("Frequency")+
    scale_x_discrete(limits = c("PriorLabels","Zeisel.Tax.Rank1","Zeisel.Tax.Rank2","Zeisel.Tax.Rank3","Zeisel.Tax.Rank4"), expand = c(.05, .05)) +
    ggtitle(paste(type,"Predictions Cross-Check",sep = " "))
  print(p5)
}
dev.off()
