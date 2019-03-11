
#Zeisel data - KI superset
load("Zeisel2018.seurat.Robj")

cells <- NULL
classes <- c("StriatDor","StriatVent")

for(i in 1:length(classes)){
  cells <- c(cells,rownames(zeisel@meta.data[which(zeisel@meta.data$Tissue == classes[i]),]))
}

zeisel.striatum <- SubsetData(object = zeisel, cells.use = cells,do.clean=T)

#Gokce 2016

load("Gokce2016.seurat.Robj")
#Saunders data
load("P60Striatum.seurat.Robj")


str <- MergeSeurat(object1 = Gokce2016, object2 = zeisel.striatum)
str <- FindVariableGenes(str, do.plot = F, display.progress = F)

hv.genes <- head(rownames(str@hvg.info), 1000)
str <- ScaleData(str, do.par=T, num.cores = 8)
str <- RunPCA(str, pc.genes = hv.genes, do.print = FALSE, pcs.compute=10)
str <- RunTSNE(str, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)


pdf(paste("Str.plots.pdf",sep=""),width=8,height = 8)
TSNEPlot(str, do.label = TRUE)
dev.off()

load("~/Documents/Harvard_Informatics/Data_Explore/mouse/FullSets/GEO/Zeisel2018/str.seurat.Robj")
TSNEPlot(str, do.label = TRUE,group.by="Tissue")
str@meta.data
