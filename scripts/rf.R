.libPaths(c(.libPaths(),"/n/home13/yasinkaymaz/biotools/Rlibs/"))
library(randomForest)

#library(flexdashboard)
library(Matrix)
library(Seurat)
library(dplyr)
#library(pheatmap)
#library(sva)
#library(d3heatmap)
#library(plotly)
#library(ggplot2)
#library(dygraphs)
#library(rbokeh)
#library(highcharter)
#library(DT)
#library(threejs)

#MODEL TRAINING
#************************************#
#load training dataset
load("~/codes/test/pbmc3k_final.Rda")
load("~/codes/test/pancreas.Robj")
setwd("~/codes/test/")
#plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
#plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = TRUE)
#plot_grid(plot1, plot2)
#ggplotly(plot1)



#Prepare training dataset:
dim(pbmc@data)
exp.pbmc <- as.data.frame(as.matrix(pbmc@data))
#It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
rownames(exp.pbmc) <- make.names(rownames(exp.pbmc))
celltype.pbmc <- as.data.frame(pbmc@meta.data[,"ClusterNames_0.6"])
rownames(celltype.pbmc) <- rownames(pbmc@meta.data)
colnames(celltype.pbmc) <- c("CellType")
head(celltype.pbmc)
celltype.pbmc <- t(celltype.pbmc)
exp.pbmc <- rbind(exp.pbmc, celltype.pbmc)
dim(exp.pbmc)

dim(pancreas@data)
exp.pancreas <- as.data.frame(as.matrix(pancreas@data))
rownames(exp.pancreas) <- make.names(rownames(exp.pancreas))
celltype.pancreas <- as.data.frame(pancreas@meta.data[,"assigned_cluster"])
rownames(celltype.pancreas) <- rownames(pancreas@meta.data)
colnames(celltype.pancreas) <- c("CellType")
head(celltype.pancreas)
celltype.pancreas <- t(celltype.pancreas)
exp.pancreas <- rbind(exp.pancreas, celltype.pancreas)


#Create Training dataset
#exp.pancreas <- t(exp.pancreas)
#exp.pbmc <- t(exp.pbmc)
trainingData <- merge(exp.pbmc, exp.pancreas, by="row.names", all=TRUE)  # merge by row names (by=0 or by="row.names")
trainingData[is.na(trainingData)] <- 0
rownames(trainingData) <- trainingData$Row.names
trainingData <- trainingData[,-c(1)]
#trainingData <- as.data.frame(t(trainingData))
#Keep the data in matrix form, otherwise randomForest will throw error: 'Error: protect(): protection stack overflow'
trainingData <- as.data.frame(t(trainingData))
#It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
#names(trainingData) <- make.names(names(trainingData))
#Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
#trainingData$CellType <- factor(trainingData$CellType)

#Optional: split the data into test and training, 20% and 80% respectively.

tsub <- trainingData[,!(names(trainingData) %in% c("CellType"))]
indx <- sapply(tsub, is.factor)
tsub[indx] <- lapply(tsub[indx], function(x) as.numeric(as.character(x)))

trainingData <- cbind(tsub, trainingData[,"CellType"])
tmat <- as.matrix(trainingData)
gc()
options("expression" =500000)

rf = randomForest(CellType~., data = tmat, norm.votes = TRUE,importance=TRUE, proximity = TRUE, ntree=100)

# TRAIN THE RF MODEL
#rf = randomForest(CellType~., data = trainingData, norm.votes = TRUE,importance=TRUE, proximity = TRUE)
print(rf)
save(rf, file="RF_model_v2")


#exp.pbmc <- t(as.data.frame(as.matrix(pbmc@data)))
#head(exp.pbmc)

#CellTypeLabels <- pbmc@meta.data[,c("orig.ident","ClusterNames_0.6")]
#traininData <- cbind(exp.pbmc, CellTypeLabels)
#It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
#names(traininData) <- make.names(names(traininData))
#Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
#traininData$ClusterNames_0.6 <- factor(traininData$ClusterNames_0.6)

#options("exppression" = 500000)
#"Rscript --max-ppsize=500000"
#print("Now training the model...")
# TRAIN THE RF MODEL
#rf = randomForest(ClusterNames_0.6~., data = traininData, norm.votes = TRUE,importance=TRUE, proximity = TRUE)
#print(rf)


#rf_importanceTable <- rf$importance
#varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")


#save(rf, file="RF_model")

#-------------------------#
#load("RF_model")

# TEST Data preparation:
# load data
#seqwell.data <- read.table(file = paste0("~/Downloads/IntegratedAnalysis_ExpressionMatrices/", 
 #                                        "pbmc_SeqWell.expressionMatrix.txt"))
#tenx.data <- read.table(file = paste0("~/Downloads/IntegratedAnalysis_ExpressionMatrices/", 
#                                      "pbmc_10X.expressionMatrix.txt"))

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

#seqwell <- CreateSeuratObject(raw.data = seqwell.data)
#seqwell <- NormalizeData(object = seqwell)
#seqwell <- ScaleData(object = seqwell)
#seqwell <- FindVariableGenes(object = seqwell, do.plot = FALSE)
#pbmc33k.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/hg19/")

# Create the object and set the initial identities
#pbmc33k  <- CreateSeuratObject(raw.data = pbmc33k.data, min.cells = 3, project = "10X_PBMC33K", names.field = 2, names.delim = "\\-")

#tenx <- CreateSeuratObject(raw.data = tenx.data)
#tenx <- NormalizeData(object = tenx)
#tenx <- ScaleData(object = tenx)
#tenx <- FindVariableGenes(object = tenx, d
