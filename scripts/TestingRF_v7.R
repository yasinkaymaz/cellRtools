.libPaths("~/biotools/Rlibs")
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")
source("main/test-functions.R")

#Testing on same data
pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
                           dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
pbmc_small <- BuildClusterTree(pbmc_small, reorder.numeric = TRUE, do.reorder = TRUE)
dim(pbmc_small@data)
tre <- pbmc_small@cluster.tree[[1]]


pbmc_small.sm <- pbmc_small
rf.sm <- CellTyperTrainer(ExpressionData = as.matrix(pbmc_small.sm@data),CellLabels = pbmc_small.sm@meta.data$tree.ident,run.name = "single.model",do.splitTest = N,PCs = 10,improve = F)


pbmc_small.sm <- CellTyper(SeuratObject = pbmc_small.sm, model = rf.sm )
head(pbmc_small.sm@meta.data)
table(pbmc_small.sm@meta.data[,c("tree.ident", "Prediction")])

pbmc_small.hm <- pbmc_small
CLoc.list <- AssessNodes(pbmc_small.hm) #Creates hierarchical random forest classifiers
#tre <- pbmc_small.hm@cluster.tree[[1]]
pbmc_small.hm <- HTyper2(SeuratObject = pbmc_small.hm, models = CLoc.list)
head(pbmc_small.hm@meta.data)
table(pbmc_small.hm@meta.data[,c("tree.ident", "Decision")])
table(pbmc_small.hm@meta.data[,c("tree.ident", "HRFPrediction")])




plot(pbmc_small.hm@meta.data$`4_classProb`, pbmc_small.sm@meta.data$`4`, pch=20,ylim = c(0,1),xlim = c(0,1))
abline(0,1)

#Testing on larger pbmc data
#pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
load("~/Documents/RFTyper/pbmc3k_final.Rda")

pbmc.sm <- pbmc
rf.sm <- CellTyperTrainer(ExpressionData = as.matrix(pbmc.sm@data),CellLabels = pbmc.sm@meta.data$tree.ident,run.name = "sm2",do.splitTest = F,PCs = 10,improve = F)

pbmc.sm <- CellTyper(SeuratObject = pbmc.sm, model = rf.sm )
head(pbmc.sm@meta.data)
table(pbmc.sm@meta.data[,c("tree.ident", "Prediction")])


pbmc.hm <- pbmc
pbmc.hm <- BuildClusterTree(pbmc.hm,do.reorder = T,reorder.numeric = T)
#PlotClusterTree(object = pbmc)
tre <- pbmc.hm@cluster.tree[[1]]
CLoc.list <- AssessNodes(pbmc.hm) #Creates hierarchical random forest classifiers
#tre <- pbmc.hm@cluster.tree[[1]]
save(CLoc.list, file="~/Documents/RFTyper/HRF_Test1/CLoc.list.pbmc3K.Rdata")


pbmc.hm <- HTyper2(SeuratObject = pbmc.hm, models = CLoc.list)

head(pbmc.hm@meta.data)
table(pbmc.hm@meta.data[,c("tree.ident", "Decision")])
table(pbmc.hm@meta.data[,c("tree.ident", "HRFPrediction")])

plot(pbmc.hm@meta.data$`8_classProb`, pbmc.sm@meta.data$`8`, pch=20,ylim = c(0,1),xlim = c(0,1))
abline(0,1)

load("~/data/Tasic2018.seurat.Robj")
Tasic2018 <- BuildClusterTree(Tasic2018,do.reorder = T,reorder.numeric = T)
PlotClusterTree(object = Tasic2018)
Tasic2018 <- HTyper2(SeuratObject = Tasic2018, models = CLoc.list)



zeisel5000subsampled <- get(load("~/codes/test/RF/HT_tests/zeisel5000subsampled.seurat.Robj"))
zeisel5000subsampled <- BuildClusterTree(zeisel5000subsampled,do.reorder = T,reorder.numeric = T)
CLoc.list <- AssessNodes(zeisel5000subsampled) #Creates hierarchical random forest classifiers
#tre <- pbmc.hm@cluster.tree[[1]]
save(CLoc.list, file="~/codes/CLoc.list.zeisel5K.Rdata")
zeisel5000subsampled <- HTyper2(SeuratObject = zeisel5000subsampled, models = CLoc.list)
table(zeisel5000subsampled@meta.data[,c("tree.ident", "Decision")])
table(zeisel5000subsampled@meta.data[,c("tree.ident", "HRFPrediction")])

rf.sm <- CellTyperTrainer(ExpressionData = as.matrix(zeisel5000subsampled@data),CellLabels = zeisel5000subsampled@meta.data$tree.ident,run.name = "sm2",do.splitTest = F,PCs = 10,improve = F)
zeisel5000subsampled <- CellTyper(SeuratObject = zeisel5000subsampled, model = rf.sm )
head(zeisel5000subsampled@meta.data)
table(zeisel5000subsampled@meta.data[,c("tree.ident", "Prediction")])

for(i in GetAllInternalNodes(tre)){print(length(CLoc.list[[as.character(i)]]$finalModel$xNames))}




source("main/test-functions.R")

ptm <- proc.time()
BuildRFClassifier <- function(object,training.genes = NULL,training.classes = NULL,verbose = TRUE,node=node,...) {
  n_tree=500
  PackageCheck('ranger')
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(
    x = training.genes,
    default = rownames(x = object@data)
  )
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = object@data[training.genes, ]
      )
    )
  )
  
  #training.data$class <- factor(x = training.classes)
  training.data <- prepareDataset(ExpressionData = as.matrix(object@data), CellLabels = training.classes, PCs = 5,run.name = node)
  
  if (verbose) {
    message("Training Classifier ...")
  }
  #classifier <- ranger::ranger(data = training.data,dependent.variable.name = "class",classification = TRUE,write.forest = TRUE,...)
  train_control <- caret::trainControl(method="cv", number=2, savePredictions = TRUE)#YK
  model.method="rf"
  names(training.data) <- make.names(names(training.data))
  #classifier <- caret::train(class~., data=training.data, trControl=train_control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)
  classifier <- caret::train(CellType~., data=training.data, trControl=train_control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=n_tree)
  return(classifier)
}

CLoc.list.n <- AssessNodes(pbmc_small.hm) #Creates hierarchical random forest classifiers
proc.time() - ptm
