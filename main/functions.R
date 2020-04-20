#Todo:
#1. implement test, validation, CV, accuracy --> DONE.
#2. Find a way to best choose KLe*Diff threshold
#3. Associate confidence with predictions. and label.
#4. Limit recursive model improvement to only one cycle in order to prevent overfitting. or instead implement a doublet detection to filter out bad cell inputs from the training data.

library(here)
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyverse)

#' A function to downsample Seurat object based on cell type identity
#' @param SeuratObj a Seurat S4 object.
#' @param IdentityCol the column 'number' in the metadata slot showing the cell type identities.
#' @usage pbmc1 <- DownSizeSeurat(SeuratObj = pbmc, IdentityCol = 7)
DownSizeSeurat <- function(SeuratObj, IdentityCol, min_n=NULL){
  cells <- NULL
  classes <- table(SeuratObj@meta.data[,IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    if( classes[type] > min_n ){
      cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = min_n, replace = F))
    }else{
      cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = classes[type], replace = F))
    }
  }
  if(class(SeuratObj)[1] == "seurat"){
    downSobj <- SubsetData(object = SeuratObj, cells.use = cells, do.clean=T)
  }else if(class(SeuratObj)[1] == "Seurat"){
    downSobj <- subset(x = SeuratObj, cells = cells)
  }
  return(downSobj)
}


#' A function to downsample a refdata table based on classLabels
#' @param RefData a data table with features as columns (last column being ClassLabels), instances in the rows.
#' @param IdentityCol the name of the column in the refdata storing class labels. default is "ClassLabels"
#' @param min_n min number of samples to downsample each class. default is the size of the minority class.
#' @usage RefData_d <- DownSampleRef(RefData = RefData)
DownSampleRef <- function(RefData, IdentityCol="ClassLabels", min_n=NULL){
  samples <- NULL
  classes <- table(RefData[, IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    if( classes[type] > min_n ){
      samples <- c(samples, sample(rownames(RefData[which(RefData[, IdentityCol] == type), ]), size = min_n, replace = F))
    }else{
      samples <- c(samples, sample(rownames(RefData[which(RefData[, IdentityCol] == type), ]), size = min_n, replace = T))
    }
  }
  RefData_d <- RefData[samples, ]
  return(RefData_d)
}




QuickSeurat <- function(SeuratObj, scale.only.var=T, PCs=20, perp=30, vars2reg) {

  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)
  if (scale.only.var == TRUE) {
    if (!missing(vars2reg)) {
      SeuratObj <- ScaleData(SeuratObj, vars.to.regress = vars2reg, genes.use = hv.genes, do.par=T, num.cores = 8)
    }else{
      SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
    }
  }else{
      if (!missing(vars2reg)) {
        SeuratObj <- ScaleData(SeuratObj, vars.to.regress = vars2reg, do.par=T, num.cores = 8)
  }else{
        SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
    }
  }
  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE, check_duplicates = FALSE)

  return(SeuratObj)

}

SeuratWrapper <- function(ExpData, ProjectLabel, NewMeta, Normalize=T, suppressLog=F, scale.only.var=T, PCs=20, perp=30, dump.files=F, min.cells=0, min.genes=0) {

  if (Normalize == TRUE) {print("Assuming the input is in count ...")
    }else{
      print("Assuming the input is in TPM ...")
      if (suppressLog == TRUE) {
        print("not taking log ...")
      }else{
        ExpData <- log1p(ExpData)
      }
    }

  SeuratObj <- CreateSeuratObject(raw.data = ExpData, project = ProjectLabel, min.cells=min.cells, min.genes = min.genes)

  if (Normalize == TRUE) {
    SeuratObj <- NormalizeData(object = SeuratObj)
    }else{
    print("Not normalizing the data since TPM is assumed ... ")
    }

  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)

  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)

  if (scale.only.var == TRUE) {
    SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
  }else{
    SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
  }

  if (!missing(NewMeta)) {
    SeuratObj <- AddMetaData(SeuratObj, NewMeta[rownames(SeuratObj@meta.data), ])
  }else{
    print("No new meta file is provided. Skipping...")
  }

  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)

  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)

  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE,check_duplicates = FALSE, perplexity=perp)

  pdf(paste(ProjectLabel,".plots.pdf", sep=""), width=8, height = 8)
  PCAPlot(SeuratObj, dim.1 = 1, dim.2 = 2)
  PCElbowPlot(SeuratObj, num.pc = PCs)
  TSNEPlot(SeuratObj, do.label = TRUE)
  dev.off()

  if (dump.files == T) {
  #Export the tSNE coordinates along with the Cluster assignment IDs
  rownames_to_column(as.data.frame(SeuratObj@dr$tsne@cell.embeddings))  %>%
    as.tibble() %>%
    add_column(Clusters=SeuratObj@meta.data$res.1) %>%
    dplyr::rename(Cellname = rowname) %>%
    as_data_frame() %>% write_csv(paste(ProjectLabel, "_tSNECoordinates_Clusters.csv", sep=""))

  #Export Normalized and Scaled Expression matrix for cells and genes in the analysis
  rownames_to_column(as.data.frame(as.matrix(SeuratObj@data))) %>%
    dplyr::rename(GeneName = rowname) %>%
    as_data_frame() %>%
    write_delim(paste(ProjectLabel, "_Normalized_Expression_matrix.txt", sep=""))
  }

  save(SeuratObj, file=paste(ProjectLabel, ".seurat.Robj", sep=""))
  return(SeuratObj)
}

SeuratCCAmerger <- function(listofObjects) {
  # Determine genes to use for CCA, must be highly variable in at least 2 datasets
  #ob.list <- list(zeisel, romanov, tasic, marques)
  ob.list <- listofObjects
  genesuse <- c()
  ids=NULL
  for (i in 1:length(ob.list)) {
    genesuse <- c(genesuse, head(rownames(ob.list[[i]]@hvg.info), 1000))
    ob.list[[i]]@meta.data$dataSource <- paste("id", i, sep="")
    ids <- c(ids, paste("id", i, sep=""))
  }
  genesuse <- names(which(table(genesuse) > 1))
  for (i in 1:length(ob.list)) {
    genesuse <- genesuse[genesuse %in% rownames(ob.list[[i]]@scale.data)]
  }

  if (length(ob.list) > 2) {
    # Run multi-set CCA
    integrated <- RunMultiCCA(ob.list, genes.use = genesuse, num.ccs = 15, add.cell.ids = ids)
    # Run rare non-overlapping filtering
    integrated <- CalcVarExpRatio(object = integrated, reduction.type = "pca", dims.use = 1:10, grouping.var = "dataSource")
    integrated <- SubsetData(integrated, subset.name = "var.ratio.pca", accept.low = 0.5)
  }else{
    #integrated <- RunCCA(object = ob.list[[1]], object2 = ob.list[[2]], genes.use = genesuse, num.cc = 15, add.cell.id = ids)
    integrated <- RunCCA(object = ob.list[[1]], object2 = ob.list[[2]], genes.use = genesuse, num.cc = 15)
  }
  # Alignment
  integrated <- AlignSubspace(integrated, reduction.type = "cca", dims.align = 1:10, grouping.var = "dataSource")
  # t-SNE and Clustering
  integrated <- FindClusters(integrated, reduction.type = "cca.aligned", dims.use = 1:10, save.SNN = T, resolution = 0.4)
  integrated <- RunTSNE(integrated, reduction.use = "cca.aligned", dims.use = 1:10)
  save(integrated, file="integrated.Aligned.seurat.Robj")
  return(integrated)
}


#A function from Hu et al.:
library("reshape2")
plot_vln <- function(t, my.genes3) {
    d <- as.matrix(t@data[intersect(my.genes3, rownames(t@data)), ])
    dd <- melt(d, id = row.names)
    dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
    dd$tree.ident <- t@ident[dd$cell]
    str(dd$tree.ident)
    dd$gene <- factor(dd$gene, levels = intersect(my.genes3, rownames(t@data)))
    ggplot(dd, aes(tree.ident, value, fill = tree.ident)) + geom_violin(scale = "width",
        trim = T, alpha = 0.8, adjust = 1) + facet_wrap(~gene, scales = "free_y",
        ncol = 1, strip.position = "right") + theme(strip.background = element_blank(),
        strip.placement = "outside", axis.text.y = element_blank(), axis.title.y = element_blank(),
        strip.text.y = element_text(colour = "red", angle = 360, size = 10),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.9)),
            legend.position = "none") + xlab("")
}
#plot_vln(scNuc2,c("Slc17a7","Gad2","Enpp6","Mog","Pdgfra","Gja1","Ctss","Flt1")) + scale_fill_manual(values=as.character(df.col$col))



library("gplots")
my.colours = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFFF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
plot_heatmap = function(t, my.genes7, my.colours = my.colours, COL = T, ROW = T,
    DEND = "none") {
    my.genes6 <- intersect(unique(my.genes7), rownames(t@data))
    Mat <- t@data[unique(my.genes6), ]
    Mat <- as.data.frame(as.matrix(Mat))
    Mat$gene <- rownames(Mat)
    Mat <- melt(Mat, id = "gene")
    Mat$cluster <- t@ident[Mat$variable]
    Mat <- Mat %>% group_by(gene, cluster) %>% dplyr::summarise(meanExp = mean(value)) %>%
        ungroup
    Mat <- as.data.frame(Mat)
    Mat <- dcast(Mat, gene ~ cluster, value.var = "meanExp")
    rownames(Mat) <- Mat$gene
    Mat <- as.matrix(Mat[, -1])
    Mat <- t(scale(t(Mat)))
    Mat <- Mat[unique(my.genes6), levels(t@ident)]
    Mat <- na.omit(Mat)
    heatmap.2(Mat, Colv = COL, Rowv = ROW, dendrogram = DEND, scale = "none",
        trace = "none", density.info = "none", col = my.colours)
}
