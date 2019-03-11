#TestingRF_v6.R
library(Seurat)

unloadNamespace("httr")
unloadNamespace("R6")
library(loomR)

lfile <- connect(filename="l5_all.loom", mode = "r+")

#guids <- paste0("gene0000", 1:length(x = lfile$row.attrs$Gene[]))#Change
guids <- make.unique(names = lfile$row.attrs$Gene[])

lfile$add.row.attribute(list(guids = guids), overwrite = TRUE)
gdf = lfile$get.attribute.df(MARGIN = 1, attribute.names = names(lfile$row.attrs), row.names = "guids")


zeisel.data  <- t(lfile$matrix[,])
rownames(zeisel.data) <- rownames(gdf)
colnames(zeisel.data) <- rownames(df)
source("~/Dropbox/codes/cellRtools/main/functions.R")
source("~/codes/cellRtools/main/functions.R")

zeisel <- SeuratWrapper(ExpData = zeisel.data, ProjectLabel = "zeisel", NewMeta = df, Normalize = T, scale.only.var = F, PCs = 100)

save(zeisel, file="Zeisel2018.seurat_updated-11-20.Robj")


load("Zeisel2018.seurat.Robj")
zeisel = SetAllIdent(zeisel, id = "ClusterName")
zeisel <- BuildClusterTree(zeisel)

pdf("Zeisel.265.clusters.tree.pdf",width = 40,height = 10)
PlotClusterTree(zeisel)
dev.off()

#node.scores <- AssessNodes(zeisel)

#Rank1

table(zeisel@meta.data$TaxonomyRank1)

cells <- NULL
classes <- names(table(zeisel@meta.data$TaxonomyRank1))
min_n <- min(table(zeisel@meta.data$TaxonomyRank1))
for(i in 1:length(classes)){
  cells <- c(cells,sample(rownames(zeisel@meta.data[which(zeisel@meta.data$TaxonomyRank1 == classes[i]),]),min_n))
}

zeisel.rank1.sub <- SubsetData(object = zeisel, cells.use = cells,do.clean=T)

zeisel.rank1.rf <- CellTyperTrainer(ExpressionData = zeisel.rank1.sub@data, CellLabels = zeisel.rank1.sub@meta.data$TaxonomyRank1, run.name = "zeisel.rank1.rf",do.splitTest = F,PCs = 100,improve = T)


#Rank2

cells <- NULL
classes <- names(table(zeisel@meta.data$TaxonomyRank2))
min_n <- min(table(zeisel@meta.data$TaxonomyRank2))
for(i in 1:length(classes)){
  cells <- c(cells,sample(rownames(zeisel@meta.data[which(zeisel@meta.data$TaxonomyRank2 == classes[i]),]),min_n))
}

zeisel.rank2.sub <- SubsetData(object = zeisel, cells.use = cells,do.clean=T)

zeisel.rank2.rf <- CellTyperTrainer(ExpressionData = zeisel.rank2.sub@data, CellLabels = zeisel.rank2.sub@meta.data$TaxonomyRank2, run.name = "zeisel.rank2.rf",do.splitTest = F,PCs = 100,improve = T)

#Rank3
cells <- NULL
classes <- names(table(zeisel@meta.data$TaxonomyRank3))
min_n <- min(table(zeisel@meta.data$TaxonomyRank3))
for(i in 1:length(classes)){
  cells <- c(cells,sample(rownames(zeisel@meta.data[which(zeisel@meta.data$TaxonomyRank3 == classes[i]),]),min_n))
}

zeisel.rank3.sub <- SubsetData(object = zeisel, cells.use = cells,do.clean=T)

zeisel.rank3.rf <- CellTyperTrainer(ExpressionData = zeisel.rank3.sub@data, CellLabels = zeisel.rank3.sub@meta.data$TaxonomyRank3, run.name = "zeisel.rank3.rf",do.splitTest = F,PCs = 100,improve = T)


#Rank4
cells <- NULL
classes <- names(table(zeisel@meta.data$TaxonomyRank4))
table_n <- table(zeisel@meta.data$TaxonomyRank4)
for(i in 1:length(classes)){
  if(as.numeric(table_n[classes[i]]) > 1000){
    print(paste("downsampling",table_n[classes[i]],classes[i],"to 1000 cells",sep=" "))
    cells <- c(cells, sample(rownames(zeisel@meta.data[which(zeisel@meta.data$TaxonomyRank4 == classes[i]),]),1000))
  }else{
    cells <- c(cells, rownames(zeisel@meta.data[which(zeisel@meta.data$TaxonomyRank4 == classes[i]),]))
  }
}
zeisel.rank4.sub <- SubsetData(object = zeisel, cells.use = cells,do.clean=T)
save(zeisel.rank4.sub, file="zeisel.rank4.sub.Rdata")
rm(zeisel)
z.rank4.rfcv <- CellTyperTrainer2(ExpressionData = zeisel.rank4.sub@data, CellLabels = zeisel.rank4.sub@meta.data$TaxonomyRank4, run.name = "zeisel.rank4.rfcv", do.splitTest = F, PCs = 100, improve = T)


#test 
source("~/Dropbox/codes/cellRtools/main/functions.R")

z.rank1.rf <- get(load("/n/holylfs/LABS/informatics/yasinkaymaz/testdata/GEO/Zeisel2018/zeisel.rank1.rf.RF_model.Robj"))

saunders <- get(load("/n/holylfs/LABS/informatics/yasinkaymaz/testdata/GEO/Saunders/Saunders_seurat.Robj"))

#Subsample saunders data
cells <- NULL
classes <- names(table(saunders@meta.data$CellType))
min_n <- min(table(saunders@meta.data$CellType))
for(i in 1:length(classes)){
  cells <- c(cells,sample(rownames(saunders@meta.data[which(saunders@meta.data$CellType == classes[i]),]),min_n))
}

saunders.sub <- SubsetData(object = saunders, cells.use = cells,do.clean=T)


saunders.sub <- CellTyper(SeuratObject = saunders.sub, model = z.rank1.rf)

saunders.sub <- FindVariableGenes(saunders.sub, do.plot = F, display.progress = F)
hv.genes <- head(rownames(saunders.sub@hvg.info), 1000)
PCs=20
saunders.sub <- ScaleData(saunders.sub, genes.use = hv.genes, do.par=T, num.cores = 8)
saunders.sub <- RunPCA(saunders.sub, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
saunders.sub <- FindClusters(saunders.sub, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
saunders.sub <- RunTSNE(saunders.sub, dims.use = 1:PCs, do.fast = TRUE,check_duplicates = FALSE)


PlotPredictions(SeuratObject = saunders.sub,model = z.rank1.rf,outputFilename = "saunders.sub.predictions.w-z.rank1.rf")

saunders.sub <- CellTyper(SeuratObject = saunders.sub, model = z.rank1.rf,priorLabels = saunders.sub@meta.data$CellType)

PlotPredictions(SeuratObject = saunders.sub,model = z.rank1.rf, outputFilename = "saunders.sub.predictions.w-z.rank1.rf")

saunders.sub@meta.data %>% group_by(Prediction, CellType, PredictionStatus) %>% tally() %>% as.data.frame()

#Saunders_P60Cerebellum_ALT
saunders.cer <- get(load("~/LabSpace/testdata/GEO/Saunders/Saunders_P60Cerebellum_ALT_seurat.Robj"))

saunders.cer <- CellTyper(SeuratObject = saunders.cer, model = z.rank1.rf, priorLabels = saunders.cer@meta.data$reason, outputFilename = "saunders.z.rank1.rf")


z.rank1.rf <- get(load("/n/holylfs/LABS/informatics/yasinkaymaz/testdata/GEO/Zeisel2018/zeisel.rank1.rf.RF_model.Robj"))
saunders.cer <- get(load("~/LabSpace/testdata/GEO/Saunders/Saunders_P60Cerebellum_ALT_seurat.Robj"))
saunders.cer <- CellTyper2(SeuratObject = saunders.cer, model = z.rank1.rf)
PlotPredictions(SeuratObject = saunders.cer, model = z.rank1.rf,outputFilename = "saunders.cer.predictions.w-z.rank1.rf")

z.rank2.rf <- get(load("/n/holylfs/LABS/informatics/yasinkaymaz/testdata/GEO/Zeisel2018/zeisel.rank2.rf.RF_model.Robj"))
saunders.cer <- get(load("~/LabSpace/testdata/GEO/Saunders/Saunders_P60Cerebellum_ALT_seurat.Robj"))
saunders.cer <- CellTyper(SeuratObject = saunders.cer, model = z.rank2.rf)
PlotPredictions(SeuratObject = saunders.cer, model = z.rank2.rf,outputFilename = "saunders.cer.predictions.w-z.rank2.rf")

z.rank3.rf <- get(load("/n/holylfs/LABS/informatics/yasinkaymaz/testdata/GEO/Zeisel2018/zeisel.rank3.rf.RF_model.Robj"))
saunders.cer <- get(load("~/LabSpace/testdata/GEO/Saunders/Saunders_P60Cerebellum_ALT_seurat.Robj"))
saunders.cer <- CellTyper(SeuratObject = saunders.cer, model = z.rank3.rf)
PlotPredictions(SeuratObject = saunders.cer, model = z.rank3.rf,outputFilename = "saunders.cer.predictions.w-z.rank3.rf")

save(zeisel.rank1.sub,file="zeisel.rank1.sub.seurat.Robj")
save(zeisel.rank2.sub,file="zeisel.rank2.sub.seurat.Robj")
save(zeisel.rank3.sub,file="zeisel.rank3.sub.seurat.Robj")
save(zeisel.rank3.sub,file="zeisel.rank3.sub.seurat.Robj")


library(caret)
load("pbmc3K.trainingData.postPCA.data")
train_control <- trainControl(method="cv", number=10)
names(getModelInfo())#To check available methods
model <- train(CellType~., data=trainingData.postPCA, trControl=train_control, method="rf")
print(model)
save(model, file="PBMC3k.rfcv.model.Rdata")

setwd("~/Documents/RFTyper/test3/")

load("../pancreas.Robj")
load("PBMC3k.rfcv.model.Rdata")
rf <- get(load("pbmc3K.RF_model.Robj"))
source("~/Dropbox/codes/cellRtools/main/functions.R")
#pancreas <- CellTyper2(SeuratObject = pancreas, model = model,priorLabels = pancreas@meta.data$assigned_cluster, outputFilename = "model.pred.pancreas")
pancreas <- HTyper(SeuratObject = pancreas, models = list(rf,rf,rf,rf), priorLabels = pancreas@meta.data$assigned_cluster, outputFilename = "model.pred.pancreas")
pancreas <- CellTyper2(SeuratObject = pancreas, model = rf,priorLabels = pancreas@meta.data$assigned_cluster, outputFilename = "model.pred.pancreas")

modelss <- list(model,model,model)
for (model in modelss){
  modelname <- deparse(substitute(model))
  print(modelname)
}


head(pancreas@meta.data)

is_alluvia_form(as.data.frame(crx), axes = 1:2, silent = TRUE)

Prior <- 'assigned_cluster'
crx <- pancreas@meta.data %>% group_by(Prior, res.0.8, Intermediate, Prediction ) %>% tally() %>% as.data.frame()

ggplot(crx,aes(y = n, axis1 = Prior , axis2=res.0.8, axis3 = Intermediate, axis4 = Prediction )) +
  geom_alluvium(aes(fill = Prediction), width = 1/12) +
  geom_stratum(width = .2, fill = "grey", color = "black") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Prior", "Clusters", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
  ggtitle("Predictions Cross-Check")


load("PCA.train.data")
load("pbmc3K.trainingData.tmp.Rdata")
rownames(trainingData)

pcadata <- data.frame(pcatrain$x, CellType = trainingData$CellType)
names(pcadata)
summary(pcatrain)


models <- c("zeisel.rank1.rf", "zeisel.rank2.rf", "zeisel.rank3.rf")
md <- c(zeisel.rank1.rf, zeisel.rank2.rf, zeisel.rank3.rf)

for (i in models){ 
  print(i)
  file <- paste("Zeisel2018/",i,".RF_model.Robj",sep="")
  assign(i,get(load(file)))
}


#with rf
setwd("~/Documents/RFTyper/presentationMaterials/")
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

hdatafiles <- c("zeisel5000subsampled.rf.predictions.w-HM.prediction-crosscheck.htable.Rdata")


zeisel.rank1.rfcv <- get(load("~/codes/test/RF/zeisel.rank1.rfcv.RF_model_notImproved.Robj"))
zeisel.rank2.rfcv <- get(load("~/codes/test/RF/zeisel.rank2.rfcv.RF_model_notImproved.Robj"))
zeisel.rank3.rfcv <- get(load("~/codes/test/RF/zeisel.rank3.rfcv.RF_model_notImproved.Robj"))
zeisel.rank4.rfcv <- get(load("~/codes/test/RF/redo_rank4/rank4.cv/zeisel.rank4.rfcv.RF_model_notImproved.Robj"))
                      

modelobjects.list <- c(zeisel.rank1.rfcv,zeisel.rank2.rfcv,zeisel.rank3.rfcv,zeisel.rank4.rfcv)
modelnames.list <- c("zeisel.rank1.rfcv","zeisel.rank2.rfcv","zeisel.rank3.rfcv","zeisel.rank4.rfcv")


model <- zeisel.rank4.rfcv
  
  #Evaluate model prediction accuracy:
  conf.mat <- model$finalModel$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>% mutate(freq = 100*value/sum(value))
  class_n = length(model$finalModel$classes)
  errorSize <- as.data.frame(cbind(model$finalModel$confusion[,"class.error"], head(colSums(model$finalModel$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- getTrainPerf(model)["TrainAccuracy"]*100
  di <- round(sqrt(class_n),digits = 0)+1
  
  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) + 
    geom_tile(color = "white")+
    coord_equal()+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types")+
    scale_x_discrete(name ="Cell Types Classes")
  
  p2 <- ggplot(errorSize) + 
    geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3", width = .5)+
    geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
    scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize),name="Class % Error rate (Dots)"))+
    labs(y="Class Size in Training (Bars)",title=paste("Model Prediction Accuracy is ",round(acc,digits = 2),"%", sep=""))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+
    scale_x_discrete(name ="Cell Type Classes")
  
pdf(paste("zeisel.rank4.rfcv","Model.evaluation.pdf",sep=""),width = 12,height = 12)

  print(p1)
  print(p2)
dev.off()
  
  
library(data.tree)

data <- read.delim("TaxonomyRank.tree.zeisel.txt",header = F)

data$pathString <- paste("ZeiselCells",data$V1, data$V2, data$V3, data$V4, sep="/")
taxa <- as.Node(data)
print(taxa)
taxaClone <- Clone(taxa)
as.data.frame(taxaClone)
library(DiagrammeR)
SetGraphStyle(taxa, rankdir = "LR")
SetEdgeStyle(taxa, arrowhead = "vee", color = "grey35", penwidth = "2px")
SetNodeStyle(taxa, style = "filled,rounded", shape = "box", fillcolor = "Yellow", 
             fontname = "helvetica", fontcolor="black", tooltip = GetDefaultTooltip,width=3)
SetNodeStyle(taxa$`Immune cells`, fillcolor = "GreenYellow", penwidth = "3px")
SetNodeStyle(taxa$Neurons, fillcolor = "Thistle", penwidth = "3px")
SetNodeStyle(taxa$Glia, fillcolor = "LightBlue", penwidth = "3px")
SetNodeStyle(taxa$`Vascular cells`, fillcolor = "Red", penwidth = "3px")
plot(taxa, direction = "descend")

taxa2 <- taxa
SetGraphStyle(taxa2, rankdir = "LR")
SetEdgeStyle(taxa2, arrowhead = "vee", color = "grey35", penwidth = "2px")
SetNodeStyle(taxa2, style = "filled,rounded", shape = "box", fillcolor = "LightBlue", 
             fontname = "helvetica", fontcolor="black", tooltip = GetDefaultTooltip,width=3)
SetNodeStyle(taxa2$`Vascular cells`, fillcolor = "LightBlue", penwidth = "3px")
SetNodeStyle(taxa2$`Immune cells`, fillcolor = "LightBlue", penwidth = "3px")
SetNodeStyle(taxa2$Neurons, fillcolor = "LightBlue", penwidth = "3px")
SetNodeStyle(taxa2$Astrocytes, fillcolor = "GreenYellow", penwidth = "3px")
plot(taxa2, direction = "descend")



print(taxa, "level")
data.frame(level = Microglia$Get('level', traversal = "ancestor"))
parent(taxa$`Vascular cells`)


pdf("taxa.tree.pdf",width = 12,height = 10)
par(mar=c(2,2,2,20))
plot(as.dendrogram(taxa), center = TRUE,horiz = T,edge.root = T)
dev.off()

crx.f$pathString <- paste("ZeiselCells",crx.f$Prior, crx.f$modelname, crx.f$modelname.1,crx.f$modelname.2,crx.f$modelname.3, sep="/")
crx.f.taxa <- as.Node(crx.f)

SetGraphStyle(crx.f.taxa, rankdir = "LR")
SetEdgeStyle(crx.f.taxa, arrowhead = "vee", color = "grey35", penwidth = 2)
SetNodeStyle(crx.f.taxa, style = "filled,rounded", shape = "box", fillcolor = "LightBlue", 
             fontname = "helvetica", fontcolor="black", tooltip = GetDefaultTooltip,width=3)
SetNodeStyle(crx.f.taxa$leaves, fillcolor = "GreenYellow", penwidth = "5px")
plot(crx.f.taxa, direction = "descend")


setwd("~/Documents/RFTyper/")
load("test3/pbmc3K.trainingData.postPCA.data")
dim(trainingData.postPCA)
library(pheatmap)
pdf("heatmap.pdf",width = 5,height = 6)
pheatmap(trainingData.postPCA[sample(1:2638,size = 10),sample(1:799,size = 10)],color=greenred(100))
dev.off()
