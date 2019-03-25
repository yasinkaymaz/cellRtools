#Testing various methods
.libPaths("~/biotools/Rlibs")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/test-functions.R")

setwd("/n/home13/yasinkaymaz/codes/test/RF/HRF_test11")
#setwd("~/Documents/RFTyper/HRF_test11/")
load("~/data/pbmc3k_final.Rda")
#source("~/Dropbox/codes/cellRtools/main/functions.R")

#create the temp pca files:
mod.rf <- CellTyperTrainer2(model.method = "rf", run.name = paste("model",sep = "."),cv_k = 10, ExpressionData = as.matrix(pbmc@data),CellLabels = pbmc@meta.data$ClusterNames_0.6,do.splitTest = F,PCs = 10,improve = F)

models <- c("rf","svmRadialWeights","kknn","lda","rda","regLogistic","wsrf","lda2","stepLDA","avNNet","nnet","pcaNNet","ORFlog")
runtimetable <- NULL
for(type in models){
  print(type)
  rmod <- paste("mod",type,sep = ".")
  rtm <- proc.time()
  assign(rmod, CellTyperTrainer2(model.method = type, run.name = paste("model",sep = "."),cv_k = 10, ExpressionData = as.matrix(pbmc@data),CellLabels = pbmc@meta.data$ClusterNames_0.6,do.splitTest = F,PCs = 10,improve = F))
  print(proc.time() - rtm)
  runtimetable <- rbind(runtimetable, type=(proc.time() - rtm))
}

rownames(runtimetable) <- models
runtimetable

results <- resamples(list(RF=mod.rf,
                          SVM=mod.svmRadialWeights,
                          KNN=mod.kknn,
                          LDA=mod.lda,
                          LDA2=mod.lda2,
                          STEPLDA=mod.stepLDA,
                          AVNET=mod.avNNet,
                          NNET=mod.nnet,
                          PCANNET=mod.pcaNNet,
                          OBRF=mod.ORFlog,
                          RDA=mod.rda,
                          RLOG=mod.regLogistic,
                          WSRF=mod.wsrf))

save(runtimetable, file="pbmc3k.runtimetable.Rdata")
save(results, file="pbmc3k.results.Rdata")
summary(results)

pdf("pbmc3k.results.plots.pdf")
bwplot(results)
dotplot(results)
dev.off()
