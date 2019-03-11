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

zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))
zeisel5000subsampled <- HTyper22(SeuratObject = zeisel5000subsampled, taxTable=TaxonomyRankTable, models = models.list, priorLabels = zeisel5000subsampled@meta.data$TaxonomyRank4, outputFilename = paste("zeisel5000subsampled",".rfcv.predictions.w-HM",sep = ""))
rm(zeisel5000subsampled)

Astrocytes <- get(load("~/LabSpace/testdata/GEO/Saunders/Astrocytes.seurat.Robj"))
Astrocytes <- HTyper22(SeuratObject = Astrocytes, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Astrocytes@meta.data$Region, outputFilename = paste("Astrocytes",".rfcv.predictions.w-HM",sep = ""))
rm(Astrocytes)

Endothelial <- get(load("~/LabSpace/testdata/GEO/Saunders/Endothelial.seurat.Robj"))
Endothelial <- HTyper22(SeuratObject = Endothelial, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Endothelial@meta.data$Region, outputFilename = paste("Endothelial",".rfcv.predictions.w-HM",sep = ""))
rm(Endothelial)

FibroblastLike <- get(load("~/LabSpace/testdata/GEO/Saunders/FibroblastLike.seurat.Robj"))
FibroblastLike <- HTyper22(SeuratObject = FibroblastLike, taxTable=TaxonomyRankTable, models = models.list, priorLabels = FibroblastLike@meta.data$Region, outputFilename = paste("FibroblastLike",".rfcv.predictions.w-HM",sep = ""))
rm(FibroblastLike)

Microglia <- get(load("~/LabSpace/testdata/GEO/Saunders/Microglia.seurat.Robj"))
Microglia <- HTyper22(SeuratObject = Microglia, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Microglia@meta.data$Region, outputFilename = paste("Microglia",".rfcv.predictions.w-HM",sep = ""))
rm(Microglia)

Mural <- get(load("~/LabSpace/testdata/GEO/Saunders/Mural.seurat.Robj"))
Mural <- HTyper22(SeuratObject = Mural, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Mural@meta.data$Region, outputFilename = paste("Mural",".rfcv.predictions.w-HM",sep = ""))
rm(Mural)

Oligodendrocytes <- get(load("~/LabSpace/testdata/GEO/Saunders/Oligodendrocytes.seurat.Robj"))
Oligodendrocytes <- HTyper22(SeuratObject = Oligodendrocytes, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Oligodendrocytes@meta.data$Region, outputFilename = paste("Oligodendrocytes",".rfcv.predictions.w-HM",sep = ""))
rm(Oligodendrocytes)

Polydendrocytes <- get(load("~/LabSpace/testdata/GEO/Saunders/Polydendrocytes.seurat.Robj"))
Polydendrocytes <- HTyper22(SeuratObject = Polydendrocytes, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Polydendrocytes@meta.data$Region, outputFilename = paste("Polydendrocytes",".rfcv.predictions.w-HM",sep = ""))
rm(Polydendrocytes)

# Testing models on real mixture datasets
#Testing on negative control
MCA_Muscle <- get(load("~/LabSpace/testdata/GEO/Han2018/MCA_Muscle.seurat.Robj"))
MCA_Muscle <- HTyper22(SeuratObject = MCA_Muscle, taxTable=TaxonomyRankTable, models = models.list, priorLabels = MCA_Muscle@meta.data$Annotation, outputFilename = paste("MCA_Muscle",".rfcv.predictions.w-HM",sep = ""))
rm(MCA_Muscle)

Gokce2016 <- get(load("~/LabSpace/testdata/GEO/Gokce2016/Gokce2016.seurat.Robj"))
Gokce2016 <- HTyper22(SeuratObject = Gokce2016, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Gokce2016@meta.data$type, outputFilename = "Gokce2016.rfcv.predictions.w-HM")
rm(Gokce2016)

Harris2018 <- get(load("~/LabSpace/testdata/GEO/Harris2018/Harris2018_corrected.seurat.Robj"))
Harris2018 <- HTyper22(SeuratObject=Harris2018, taxTable=TaxonomyRankTable, priorLabels=Harris2018@meta.data$Comments, models=models.list, outputFilename="Harris2018.rfcv.predictions.w-HM")
rm(Harris2018)

Romanov2017 <- get(load("~/LabSpace/testdata/GEO/Romanov2017/Romanov2017.seurat.Robj"))
Romanov2017 <- HTyper22(SeuratObject=Romanov2017, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Romanov2017@meta.data$level1.class, outputFilename="Romanov2017.rfcv.predictions.w-HM")
rm(Romanov2017)

Hrvatin <- get(load("~/LabSpace/testdata/GEO/Hrvatin2018/Hrvatin.seurat.Robj"))
Hrvatin <- HTyper22(SeuratObject=Hrvatin, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Hrvatin@meta.data$celltype, outputFilename="Hrvatin.rfcv.predictions.w-HM")
rm(Hrvatin)

Hochgerner2018 <- get(load("~/LabSpace/testdata/GEO/Hochgerner2018/Hochgerner2018.seurat.Robj"))
Hochgerner2018 <- HTyper22(SeuratObject=Hochgerner2018, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Hochgerner2018@meta.data$res.1, outputFilename="Hochgerner2018.rfcv.predictions.w-HM")
rm(Hochgerner2018)

Campbell2017  <- get(load("~/LabSpace/testdata/GEO/Campbell2017/Campbell2017.seurat.Robj"))
Campbell2017 <- HTyper22(SeuratObject = Campbell2017, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Campbell2017@meta.data$res.1, outputFilename = "Cambell2017.rfcv.predictions.w-HM")
rm(Campbell2017)

Tasic2016  <- get(load("~/LabSpace/testdata/GEO/Tasic2016/Tasic2016.seurat.Robj"))
Tasic2016 <- HTyper22(SeuratObject=Tasic2016, taxTable=TaxonomyRankTable, models = models.list, priorLabels = Tasic2016@meta.data$res.1,  outputFilename="Tasic2016.rfcv.predictions.w-HM")
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
  crx.f <- crx %>% mutate(freq = n*100 / sum(n)) %>% filter(freq > 1)
  p5 <- ggplot(crx.f,aes_string(y = "n", axis1 = names(crx.f)[1], axis2 = names(crx.f)[2], axis3 = names(crx.f)[3], axis4 = names(crx.f)[4], axis5 = names(crx.f)[5], axis6 = names(crx.f)[6])) +
    geom_alluvium(aes_string(fill = names(crx.f)[6]), width = 0, knot.pos = 1/4) +
    guides(fill = FALSE)+
    geom_stratum(width = 1/12, fill = "grey", color = "red") +
    geom_label(stat = "stratum", label.strata = TRUE ) +
    ylab("Frequency")+
    scale_x_discrete(limits = c("PriorLabels","Zeisel.Tax.Rank1","Zeisel.Tax.Rank2","Zeisel.Tax.Rank3","Zeisel.Tax.Rank4","FinalPrediction"), expand = c(.05, .05)) +
    ggtitle(paste(type,"Predictions Cross-Check",sep = " "))
  print(p5)
}
dev.off()



freq.Final <- zeisel5000subsampled@meta.data[,c("TaxonomyRank4","FinalPrediction")] %>% as.tibble() %>% dplyr::count(TaxonomyRank4, FinalPrediction) %>% group_by(TaxonomyRank4) %>% mutate(prop = 100*prop.table(n)) %>% as.data.frame()

pp <- ggplot(freq.Final, aes(TaxonomyRank4, FinalPrediction, fill=prop)) +
    geom_tile(color = "white")+
    coord_equal()+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types with Final Model (HRF)")+
    scale_x_discrete(name ="Cell Types Classes")

save_plot(filename = paste("ZeiselTax4-FinalPrediction",".CrossCheck.pdf",sep=""),plot = pp,base_height = 12, base_width = 12)
freq.Final %>% write_tsv(path="ZeiselTax4-FinalPrediction.CrossCheck.txt")


fullHtable <- get(load("zeisel5000subsampled.rfcv.predictions.w-HM.prediction-crosscheck.Full_htable.Rdata"))

freq.mod3 <- fullHtable %>% as.tibble() %>% dplyr::count(Prior, modelname.3) %>% group_by(Prior) %>% mutate(prop = 100*prop.table(n)) %>% as.data.frame()

pp <- ggplot(freq.mod3, aes(Prior, modelname.3, fill=prop)) +
    geom_tile(color = "white")+
    coord_equal()+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types with model-4")+
    scale_x_discrete(name ="Cell Types Classes")

save_plot(filename = paste("ZeiselTax4-Model4",".-CrossCheck.pdf",sep=""),plot = pp,base_height = 12, base_width = 12)
freq.mod3 %>% write_tsv(path="ZeiselTax4-Model4.CrossCheck.txt")
