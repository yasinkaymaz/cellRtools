#Todo:
#1. implement test, validation, CV, accuracy --> DONE.
#2. Find a way to best choose KLe*Diff threshold
#3. Associate confidence with predictions. and label.
#4. Limit recursive model improvement to only one cycle in order to prevent overfitting. or instead implement a doublet detection to filter out bad cell inputs from the training data.

library(here)

RunGSEAforClusters <- function(SeuratObj, Cond1, Cond2, GeneSet=here('data/GeneSetDatabases/MousePath_GO_gmt.gmt'), outputDir=getwd(), ...){

  SeuratObj = SetAllIdent(SeuratObj, id = 'cluster.states')
  clusters = sort(unique(as.numeric(SeuratObj@meta.data$tree.ident)))
  pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)

  for (i in clusters) {
    #Create GSEA input files
    filePrefix= paste(Cond1,Cond2,"Cluster",i,"ExpMatrix",sep="_")
    CreateGSEAinput(SeuratObj = SeuratObj, Cond1 = Cond1, Cond2 = Cond2, clusterid = i, filePrefix = filePrefix, ...)
    #Run GSEA for the cluster i comparing cells from Cond1 and Cond2
    RunGSEA(InputPrefix = filePrefix, GeneSet=GeneSet, ...)

    setTxtProgressBar(pb,i)
    print(c('GSEA with Cluster', i,' is completed.'))
  }
}

RunDGEA <- function(SeuratObj, Cond1, Cond2, outputDir=getwd(), ...){
  #First install MAST:
  #install.packages("BiocManager")
  #BiocManager::install("MAST")
  #library(scater)
  library(dplyr)
  SeuratObj = SetAllIdent(SeuratObj, id = 'cluster.states')
  DEGs = list()
  clusters = sort(unique(as.numeric(SeuratObj@meta.data$tree.ident)))
  pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
  for (i in clusters) {
    ident.1 = paste(Cond1,"-",i, sep = '')
    ident.2 = paste(Cond2,'-',i, sep = '')
    a =FindMarkers(SeuratObj, ident.1 = ident.1, ident.2 = ident.2, test.use = 'MAST', only.pos = FALSE, latent.vars = c('nUMI', 'percent.mito'), logfc.threshold = 0.1, min.pct = 0.05)
    a$gene = rownames(a)
    a$cluster = paste(i,sep = '')
    DEGs = bind_rows(DEGs,a)
    setTxtProgressBar(pb,i)
    print(c('Finished with Cluster', i))
  }
  save(DEGs, file = paste(outputDir,"/",Cond1,"-",Cond2,"DEGs.Rdata",sep = ""))
  #unloadNamespace(ns = 'scater')
  return(DEGs)
}

CreateGSEAinput <- function(SeuratObj, Cond1, Cond2, outputDir=getwd(), clusterid, filePrefix, ...){

  ident.1 = paste(Cond1,"-",clusterid, sep = '')
  ident.2 = paste(Cond2,'-',clusterid, sep = '')
  #Subset expression data into cells from two conditions:
  sub.data.id1 <- as.data.frame(as.matrix(x = SeuratObj@data[, WhichCells(object = SeuratObj, ident = ident.1)]))
  sub.data.id2 <- as.data.frame(as.matrix(x = SeuratObj@data[, WhichCells(object = SeuratObj, ident = ident.2)]))
  #Store cell numbers in each condition here:
  c1 <- dim(sub.data.id1)[2]
  c2 <- dim(sub.data.id2)[2]
  #merge the datasets from both conditions into a df
  tpm <- cbind( sub.data.id1,sub.data.id2)
  new.df <- cbind(Description="",tpm)
  #get the matrix dimensions
  dimensions <- dim(tpm)

  #Create a GCT file: for details: https://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT
  header1="#1.2"
  header2=paste(dimensions[1],dimensions[2],sep="\t")
  write.table(header1, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE,col.names=FALSE,row.names=FALSE   )
  write.table(header2, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE  )
  write.table(new.df, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE, append=TRUE,col.names=NA   )

  #create a CLS file: for details: https://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS
  conditions = c(ident.1,ident.2)
  header=paste(dimensions[2], "2", "1",sep=" ")
  line2=paste("#",conditions[1],conditions[2], sep=" ")
  line3=paste( rep(c(conditions[1],conditions[2]), c(c1,c2)),sep = " " )
  write.table(header, file=paste(filePrefix,".cls",sep=""), sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
  write.table(line2,file=paste(filePrefix,".cls",sep=""), sep=" ", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE)
  linex=line3[1]
  for (i in 2:length(line3)){
    linex <- paste(linex,line3[i],sep =" ")}
  write.table(linex,file=paste(filePrefix,".cls",sep=""), sep=" ", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE)
}

RunGSEA <- function(InputPrefix, GeneSet, outputDir=getwd(), ...){
  GSEA.program.location <- paste(srcdir,"/GSEA.1.1.R",sep="")
  source(GSEA.program.location, verbose=F, max.deparse.length=9999)

  doc.STRING= paste(InputPrefix,substr(GeneSet,1,nchar(GeneSet)-4), sep="_")
  print(doc.STRING)
  print(InputPrefix)

  GSEA(   # Input/Output Files :-------------------------------------------
          input.ds =  paste(outputDir,"/",InputPrefix,".gct",sep = ""),           # Input gene expression Affy dataset file in RES or GCT format
          input.cls = paste(outputDir,"/",InputPrefix,".cls",sep = ""),           # Input class vector (phenotype) file in CLS format
          gs.db =   paste(srcdir,"/../GeneSetDatabases/",GeneSet,sep=""),         # Gene set database in GMT format
          output.directory      = paste(outputDir,"/",sep = ""),        # Directory where to store output and results (default: "")
          #  Program parameters :-------------------------------------------------------------------------------------------------------------------------
          doc.string            = doc.STRING,   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
          non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
          reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
          nperm                 = 1000,            # Number of random permutations (default: 1000)
          weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
          nom.p.val.threshold   = 0.001,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
          fwer.p.val.threshold  = 0.001,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
          fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
          topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
          adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
          gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
          gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
          reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
          preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
          random.seed           = 3338,            # Random number generator seed. (default: 123456)
          perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
          fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
          replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
          save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
          OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
          use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
  )

}

SummarizeGSEAoutputs <- function(GSEAoutputDir="./"){
  #This function returns the table of all Enrichment results with corrected p-values.
  library(tidyverse)
  setwd(GSEAoutputDir)
  
  majorSummaryTable <- NULL
  GSreportsTable <- NULL
  mySumFiles <- list.files(pattern="*SUMMARY.RESULTS.REPORT*")

  for (i in 1:length(mySumFiles)){

    sumTable <- read.delim(mySumFiles[i]) %>% as.tibble() %>% add_column(Comparison=strsplit(mySumFiles[i],"_Clust")[[1]][1],EnrichmentDirection_ClusterID=strsplit(mySumFiles[i],"\\.")[[1]][5])
    majorSummaryTable <- bind_rows(majorSummaryTable, sumTable)
    
    #for each Gene set, j, in the summary table:
    for(j in 1:length(read.delim(mySumFiles[i])[,1])){
      #the Gene Set j from the directory: Get the file prefix from the Summary file name + combine with gene set name + add ".txt" to the end.
      geneSetReportfile=list.files(pattern=paste(strsplit(mySumFiles[i],"\\.")[[1]][1], (read.delim(mySumFiles[i]) %>% as.tibble() %>% select(GS) %>% c())[[1]][j],"report",strsplit(mySumFiles[i],"\\.")[[1]][5], "*.txt", sep = "."))
      
      #if (!identical(geneSetReportfile, character(0))){
      if (!identical(geneSetReportfile, character(0)) && (geneSetReportfile != "Subordinate_Control_Cluster_19_ExpMatrix_Calvin_manual_genesets.neuromuscular junction.report.Control-19.12.txt")){
          
        gs.reporttable <-  read.delim(geneSetReportfile) %>% 
          as.tibble() %>%
          dplyr::filter(CORE_ENRICHMENT == "YES") %>% # filter out genes which are not in the Leading Edge.
          add_column(
              Comparison = strsplit(mySumFiles[i],"_Clust")[[1]][1], #Create a column for Comparison type, ex; 'Dominant_Control'
              EnrichmentDirection_ClusterID = strsplit(mySumFiles[i],"\\.")[[1]][5], #Create a column for Enrichment direction, ex; 'Control-1'. This also shows the cluster id.
              GS = (read.delim(mySumFiles[i]) %>% as.tibble() %>% select(GS) %>% c())[[1]][j] #Create a column for Gene Set name.
            )
        GSreportsTable <- bind_rows(GSreportsTable, gs.reporttable)
        
      }else{break}#closes ifelse for report file existance.
    }#closes loop for j
  }#closes loop for i
  
  majorSummaryTable <- majorSummaryTable %>% as.tibble() %>% mutate(pAdj.Nom=p.adjust(NOM.p.val,method="BH")) %>% arrange(pAdj.Nom)
  sigtable <- majorSummaryTable %>% dplyr::filter(pAdj.Nom < 0.05) %>% unite(plotfilename, Comparison,GS,EnrichmentDirection_ClusterID,sep="*",remove = FALSE)
  #Write the main table and only significant enrichments to separate files:
  majorSummaryTable %>% write_tsv(file.path("All.Enrichment.stats.txt"))
  sigtable %>% write_tsv(file.path("significant.Enrichments.txt"))
  
  sig.GSreportsTable=NULL
  for(g in 1:dim(sigtable)[1]){
    sig.g <- GSreportsTable %>% filter(Comparison == sigtable[g,]$Comparison, EnrichmentDirection_ClusterID == sigtable[g,]$EnrichmentDirection_ClusterID, GS == sigtable[g,]$GS) %>% select(-SYMBOL, -DESC) %>% separate(EnrichmentDirection_ClusterID, into=c("EnrichmentDirection","cluster"))
    sig.GSreportsTable <- bind_rows(sig.GSreportsTable, sig.g)
  }
  
  #return(majorSummaryTable)
  return(sig.GSreportsTable)
}

SeuratWrapper <- function(ExpData, ProjectLabel, NewMeta, Normalize=T, scale.only.var=T, PCs=20, perp=30, dump.files=F){
  
  if(Normalize == TRUE){print("Assuming the input is in count ...")
    }else{
      print("Assuming the input is in TPM ...")
    ExpData <- log1p(ExpData)
    }
  
  SeuratObj <- CreateSeuratObject(raw.data = ExpData, project = ProjectLabel, min.genes = 500)
  
  if (Normalize == TRUE) {
    SeuratObj <- NormalizeData(object = SeuratObj)
    }else{
    print("Not normalizing the data since TPM is assumed ... ")
    }
  
  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
  
  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)
  
  if (scale.only.var == TRUE){
    SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
  }else{
    SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
  }
  
  if(!missing(NewMeta)){
    SeuratObj <- AddMetaData(SeuratObj, NewMeta)    
  }else{
    print("No new meta file is provided. Skipping...")
  }
  
  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
  
  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
  
  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE,check_duplicates = FALSE, perplexity=perp)
  
  pdf(paste(ProjectLabel,".plots.pdf",sep=""),width=8,height = 8)
  PCAPlot(SeuratObj, dim.1 = 1, dim.2 = 2)
  PCElbowPlot(SeuratObj, num.pc = PCs)
  TSNEPlot(SeuratObj, do.label = TRUE)
  dev.off()
  
  if(dump.files == T){
  #Export the tSNE coordinates along with the Cluster assignment IDs
  rownames_to_column(as.data.frame(SeuratObj@dr$tsne@cell.embeddings))  %>%
    as.tibble() %>% 
    add_column(Clusters=SeuratObj@meta.data$res.1) %>%
    dplyr::rename(Cellname = rowname) %>%
    as_data_frame() %>% write_csv(paste(ProjectLabel,"_tSNECoordinates_Clusters.csv",sep=""))
  
  #Export Normalized and Scaled Expression matrix for cells and genes in the analysis
  rownames_to_column(as.data.frame(as.matrix(SeuratObj@data))) %>% 
    dplyr::rename(GeneName = rowname) %>% 
    as_data_frame() %>% 
    write_delim(paste(ProjectLabel,"_Normalized_Expression_matrix.txt",sep=""))
  }
  
  save(SeuratObj, file=paste(ProjectLabel,".seurat.Robj",sep=""))
  return(SeuratObj)
}

SeuratCCAmerger <- function(listofObjects){
  # Determine genes to use for CCA, must be highly variable in at least 2 datasets
  #ob.list <- list(zeisel, romanov, tasic, marques)
  ob.list <- listofObjects
  genesuse <- c()
  ids=NULL
  for (i in 1:length(ob.list)) {
    genesuse <- c(genesuse, head(rownames(ob.list[[i]]@hvg.info), 1000))
    ob.list[[i]]@meta.data$dataSource <- paste("id",i,sep="")
    ids <- c(ids, paste("id",i,sep=""))
  }
  genesuse <- names(which(table(genesuse) > 1))
  for (i in 1:length(ob.list)) {
    genesuse <- genesuse[genesuse %in% rownames(ob.list[[i]]@scale.data)]
  }
  
  if(length(ob.list) > 2){
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

FeaturePrep <- function(SeuratObj, gene.set, ScoreName){
  library(matrixStats)
  # Get mean expression of genes of interest per cell
  #mean.exp <- colMeans(x = SeuratObj@data[gene.set, ], na.rm = TRUE)
  mean.exp <- colMaxs(x = 2^as.matrix(SeuratObj@data[gene.set, ]), na.rm = TRUE)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = SeuratObj@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
        "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
    SeuratObj@meta.data[[ScoreName]] <- mean.exp
  }
  return(SeuratObj)
}

selectGenes.best.loadings <- function(trainingExpData, pcs, num, caption="Highest loading"){
  #Second version:
  #p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom
  #Weighted gene picking depending on PC number: Initial PCs give more genes.
  #For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
  print("Performing PCA...")
  #Do PCA here on input data
  library(matrixStats)
  trainingExpData <- trainingExpData[, apply(trainingExpData, 2, var) != 0]
  trainingExpData <- trainingExpData[, which(colVars(as.matrix(trainingExpData)) > 0.05)]
  
  pcatrain <- prcomp(trainingExpData, center = FALSE, scale=FALSE, rank. = pcs)
  save(pcatrain,file="PCA.train.data")
  
  print("Selecting the genes as best features...")
  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  pb <- txtProgressBar(min = 0, max = pcs, style = 3)
  for(i in 1:pcs){
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[,i]),decreasing = TRUE),i]
    Gn <- round((pcs-i+1)*(num*2)/(pcs*(pcs+1)))
    TotalGenes = as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai,Gn)),bestgenes=head(orderedpcai,Gn))
    load <- rbind(load, top)
    setTxtProgressBar(pb,i)
    cat("\n",'Picking the best genes from first', i,'PCs is completed.',"\n")
  }
  
  bestgenes <- unique(load$genes)
  return(bestgenes)
}

prepareDataset <- function(ExpressionData, CellLabels, run.name, PCs, featureGeneSet){
  #Required: This will take a Normalized expression data matrix, rows as genes and columns as cells. Example: as.matrix(SeuratObject@data)
  #Required: A list of cell labels. same dimension as colnames(input expression). Example: SeuratObject@meta.data$res.1
  
  #Transpose the matrix cols <--> rows t()
  #Keep the data in matrix form, otherwise randomForest will throw error: 'Error: protect(): protection stack overflow'
  trainingData <- as.data.frame(t(ExpressionData))
  #It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
  names(trainingData) <- make.names(names(trainingData))
  trainingData$CellType <- CellLabels
  #Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
  trainingData$CellType <- factor(trainingData$CellType)
  
  train <- droplevels(trainingData[,!(names(trainingData) %in% c("CellType"))])
  indx <- sapply(train, is.factor)
  train[indx] <- lapply(train[indx], function(x) as.numeric(as.character(x)))
 
  save(trainingData, file=paste(run.name,".trainingData.tmp.Rdata",sep = ""))
  rm(trainingData, tsub)
  gc()
  
  if(missing(featureGeneSet)){
    #Perform PCA on data (optional: only on training portion) and select genes for features
    genes <- as.character(selectGenes.best.loadings(trainingExpData = train, pcs = PCs, num = 2000))
  }else{
    genes <- make.names(featureGeneSet)
  }#closes.if.missing.featureset
  
  load(paste(run.name,".trainingData.tmp.Rdata",sep = ""))
  
  pdf(paste(run.name,"_PCAplots.pdf",sep=""),width = 10,height = 10)
  for (i in 1:PCs){
    pci <- paste("PC",i,sep="")
    pcj <- paste("PC",i+1,sep="")
    print(ggplot(trainingData, aes_string(x=pci, y=pcj, color="CellType"))+geom_point() )
    par()}
  dev.off()

  trainingData.postPCA <- droplevels(trainingData[,c(genes,"CellType")])
  
  save(trainingData.postPCA, file=paste(run.name, ".trainingData.postPCA.data", sep = ""))
  file.remove(paste(run.name, ".trainingData.tmp.Rdata", sep = ""))
  
  return(trainingData.postPCA)
}

CellTyperTrainer2 <- function(ExpressionData, CellLabels, model.method="rf", run.name, do.splitTest=F, PCs, improve.rf=F, cv_k=5){
  library(randomForest)
  library(rfUtilities)
  library(tidyverse)
  library(caret)
  
  if(missing(PCs)){
    PCs=length(unique(CellLabels))
  }else{
    PCs=PCs
  }
  
  ExpressionData <- as.matrix(ExpressionData)
  
  if(file.exists(paste(run.name,".trainingData.postPCA.data",sep = ""))){
    print("Training data already exists...")
    trainingData <- get(load(paste(run.name,".trainingData.postPCA.data",sep = "")))
  }else{
    print("creating the training data...")
    trainingData <- prepareDataset(ExpressionData = ExpressionData, CellLabels = CellLabels, PCs = PCs, run.name = run.name)
  }
  
  #k-fold Cross Validation
  train_control <- trainControl(method="cv", number=cv_k, savePredictions = TRUE)
  model <- train(CellType~., data=trainingData, trControl=train_control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)
  
  if((improve.rf == T) & (model.method = "rf")){
    
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    
    rfvotes <- as.data.frame(model$finalModel$votes)
    rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
    rfvotes$label <- apply(model$finalModel$votes, 1, function(x)  names(x)[which.max(x)] )
    rfvotes$inputnames <- rownames(rfvotes)
    e <- as.data.frame(table(rfvotes$label))
    Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
    
    #Defaults
    filter_p = 0.05
    bestvote_cutoff = 0.7
    badinput_ids <- NULL
    round_n=1
    badinput_stats <- data.frame()
    currentscore = Bestscore$median
    toss_n = dim(rfvotes)[1]
    
    #While Median Best votes score overall is larger in the new iteration AND the number of inputs need to be tossed is larger %1 of all input data, then continue.
    #while (Bestscore$median >= currentscore && toss_n > round(0.01*dim(rfvotes)[1]) ){
    while (round_n < 2 ){#run this only once...
      
      print(paste("Round number ",round_n))
      print(paste("Current score is", currentscore,". toss_n is", toss_n, ". Fractions is", round(0.01*dim(rfvotes)[1])))
      
      for(i in 1:length(e$Var1)){
        badinputs <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% select(inputnames) %>% c()
        badinputs.m <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% summarize(mean=mean(bestvote)) %>% c()
        badinput_ids <- unique(c(badinput_ids, badinputs$inputnames))
        classBestscore <- rfvotes %>% as.tibble() %>%  filter(., label == e$Var1[i]) %>% summarize(median=median(bestvote)) %>% c()
        class_badinput_stats <- data.frame(class=e$Var1[i], classInputSize=e$Freq[i], allBestscoreMedian=Bestscore$median, classBestscoreMedian=classBestscore$median, tossedInput=length(badinputs$inputnames), tossedBestvoteMean=badinputs.m$mean, iteration=round_n)
        badinput_stats <- rbind(badinput_stats, class_badinput_stats)
      }
      
      badinput_stats[is.nan(badinput_stats)] <- 0
      toss_n <- badinput_stats %>% as.tibble() %>% filter(., iteration == round_n) %>% summarise(n=sum(tossedInput)) %>% c()
      toss_n <- toss_n$n
      
      print(badinput_stats)
      
      #filter input using the bad input list generated in the previous iteration
      trainingData <- trainingData[which(!rownames(trainingData) %in% badinput_ids ),]
      
      #run the RF again with the updated training set
      model <- train(CellType~., data=trainingData, trControl=train_control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)
      
      #Check again to see if there is room to improve:
      rfvotes <- as.data.frame(model$finalModel$votes)
      rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
      rfvotes$label <- apply(model$finalModel$votes, 1, function(x)  names(x)[which.max(x)] )
      rfvotes$inputnames <- rownames(rfvotes)
      e <- as.data.frame(table(rfvotes$label))
      Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
      #update the round
      round_n = round_n + 1
    }#closes the while loop
    
  }#closes the improve option
  
  pdf(paste(run.name,".pdf",sep=""),width= 1.5*class_n, height = 1.5*class_n)
  require(gridExtra)
  errorSize <- as.data.frame(cbind(model$finalModel$confusion[,"class.error"],
                                   head(colSums(model$finalModel$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- getTrainPerf(model)["TrainAccuracy"]*100# This is the mean accuracy of each fold at model$resample
  
  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types")+
    scale_x_discrete(name ="Cell Types Classes")
  
  p2 <- ggplot(errorSize)  + 
    geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3")+
    geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
    scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize),name="Class % Error rate (Dots)"))+
    labs(y="Class Size in Training (Bars)",title=paste("Model Prediction Accuracy is ",round(acc,digits = 2),"%", sep=""))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+
    scale_x_discrete(name ="Cell Type Classes")
  
  grid.arrange(p1, p2, nrow=2)
  
  dev.off()
  
  save(model, file=paste(run.name,".RF_model.Robj",sep = ""))
  return(model)
}

PlotPredictions2 <- function(SeuratObject, model, priorLabels, outputFilename="plotpredictions"){
  #Evaluate model prediction accuracy:
  conf.mat <- model$finalModel$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>% mutate(freq = 100*value/sum(value))
  class_n = length(model$finalModel$classes)
  errorSize <- as.data.frame(cbind(model$finalModel$confusion[,"class.error"], head(colSums(model$finalModel$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- getTrainPerf(model)["TrainAccuracy"]*100
  di <- round(sqrt(class_n),digits = 0)+1
  
  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types")+
    scale_x_discrete(name ="Cell Types Classes")
  
  p2 <- ggplot(errorSize)  + 
    geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3")+
    geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
    scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize),name="Class % Error rate (Dots)"))+
    labs(y="Class Size in Training (Bars)",title=paste("Model Prediction Accuracy is ",round(acc,digits = 2),"%", sep=""))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+
    scale_x_discrete(name ="Cell Type Classes")
  
  #Prediction Performance
  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=Prediction,fill=Prediction),stat = "count")+
    geom_violin(aes(x=Prediction,y=BestVotesPercent*max(table(SeuratObject@meta.data$Prediction)),fill=Prediction))+ 
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$Prediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))
  
  p4 <- TSNEPlot(SeuratObject, group.by="Prediction",do.label=T, do.return = T)
  
  
  ff <- FeaturePlot(object = SeuratObject,
                    features.plot = model$finalModel$classes,
                    cols.use = c("grey", "blue"),
                    reduction.use = "tsne",do.return = T)
  
  q1 <- FeaturePlot(SeuratObject, features.plot = "BestVotesPercent", no.legend = F, cols.use = c("gray","red"),do.return = T)
  q2 <- FeaturePlot(SeuratObject, features.plot = "KLe", no.legend = F, cols.use = c("gray","purple"),do.return = T)
  q3 <- ggplot(SeuratObject@meta.data, aes(KLe, Diff, color= PredictionStatus))+ geom_point(size=.6)
  
  p1p2 <- cowplot::plot_grid(p1,NULL,NULL, p2, ncol = 2, nrow=2)
  save_plot(filename = paste(outputFilename,".mode.-stats.pdf",sep=""),plot = p1p2,base_height = 12, base_width = 12)
  
  p3p4 <- cowplot::plot_grid(p3,NULL,NULL,p4, ncol = 2, nrow=2)
  save_plot(filename = paste(outputFilename,".prediction-stats.pdf",sep=""),plot = p3p4,base_height = 12, base_width = 12)
  
  ff <- cowplot::plot_grid(plotlist=ff, nrow = di, ncol = di-1)
  save_plot(filename = paste(outputFilename,".prediction-probability-projection.pdf",sep=""),plot = ff,base_height = 12, base_width = 12)
  
  q1q2q3 <- cowplot::plot_grid(plotlist = list(q1$BestVotesPercent, q2$KLe, q3 ), ncol = 2, nrow = 2)
  save_plot(filename = paste(outputFilename,".prediction-quality.pdf",sep=""),plot = q1q2q3,base_height = 12, base_width = 12)
  
}

CellTyper2 <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions"){
  
  library(caret)
  library(randomForest)
  library(tidyverse)
  
  if(!missing(SeuratObject)){
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)
  colnames(testExpSet) <- make.names(colnames(testExpSet))
  #Prepare Test Expression set
  testsub <- testExpSet[,which(colnames(testExpSet) %in% model$finalModel$xNames)]
  missingGenes <- model$finalModel$xNames[which(!model$finalModel$xNames %in% colnames(testExpSet))]
  print(model$finalModel$importance[missingGenes,])
  
  missingGenes.df <- data.frame(matrix(0, ncol = length(missingGenes), nrow = length(rownames(testExpSet))))
  colnames(missingGenes.df) <- missingGenes
  TestData <- cbind(testsub, missingGenes.df)
  TestData <- TestData[,model$finalModel$xNames]
  mmDGm <- mean(model$finalModel$importance[missingGenes,])
  mmDGf <- mean(model$finalModel$importance[which(!model$finalModel$xNames %in% missingGenes),])
  
  cat("Number of Features (genes) to be considered is", length(colnames(testsub)), '\n', 
      "Number of missing Features set to zero is", length(missingGenes), '\n', 
      "Mean MeanDecreaseGini of missing genes is", mmDGm, '\n',
      "Mean MeanDecreaseGini of Featured genes is", mmDGf, '\n',
      sep = ' ')
  
  if((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )){
    warning("A significant portion of features are missing...")
    }
  
  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class_n <- length(model$finalModel$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class_n, class_n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="raw")
  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class_n) 
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class_n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class_n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class_n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()
  
  if(missing(priorLabels)){
    print("Prior class labels are not provided!")
    
  }else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    testPred <- cbind(testPred, priorLabels)
    
    #Plot the crosscheck here:
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    crx <- testPred %>% group_by(Prior, res.1, Intermediate, Prediction) %>% tally() %>% as.data.frame()
    
    p5 <- ggplot(crx,aes(y = n, axis1 = Prior , axis2 = Prediction )) +
      geom_alluvium(aes(fill = Prediction), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Prior", "Clusters", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")
    
    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class_n, base_width = 1.2*class_n)
  }
  
  if(!missing(SeuratObject)){
    #update predictions in the meta.data slot 
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)
    
    PlotPredictions2(SeuratObject = SeuratObject, model = model, outputFilename = outputFilename)
    
    return(SeuratObject)
    
  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function


####################################Previous version functions###############################################
CellTyperTrainer <- function(ExpressionData, CellLabels, run.name, do.splitTest=F, PCs, improve=T){
  library(randomForest)
  library(rfUtilities)
  library(tidyverse)
  
  if(missing(PCs)){
    PCs=length(unique(CellLabels))
  }else{
    PCs=PCs
  }
  
  ExpressionData <- as.matrix(ExpressionData)
  
  if(file.exists(paste(run.name,".trainingData.postPCA.data",sep = ""))){
    print("Training data already exists...")
    trainingData <- get(load(paste(run.name,".trainingData.postPCA.data",sep = "")))
  }else{
    print("creating the training data...")
    trainingData <- prepareDataset(ExpressionData = ExpressionData, CellLabels = CellLabels, do.splitTest = do.splitTest, PCs = PCs, run.name = run.name)
  }
  
  #Added: "sampsize=c(table(trainingData$CellType))". Revisit this later to make sure it is working as expected...
  rf <- randomForest(CellType~., data = trainingData, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData$CellType)))
  print(rf)
  
  if(improve == T){
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    
    rfvotes <- as.data.frame(rf$votes)
    rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
    rfvotes$label <- apply(rf$votes, 1, function(x)  names(x)[which.max(x)] )
    rfvotes$inputnames <- rownames(rfvotes)
    e <- as.data.frame(table(rfvotes$label))
    Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
    
    #Defaults
    filter_p = 0.05
    bestvote_cutoff = 0.7
    badinput_ids <- NULL
    round_n=1
    badinput_stats <- data.frame()
    currentscore = Bestscore$median
    toss_n = dim(rfvotes)[1]
    
    #While Median Best votes score overall is larger in the new iteration AND the number of inputs need to be tossed is larger %1 of all input data, then continue.
    while (Bestscore$median >= currentscore && toss_n > round(0.01*dim(rfvotes)[1]) ){
      
      print(paste("Round number ",round_n))
      print(paste("Current score is", currentscore,". toss_n is", toss_n, ". Fractions is", round(0.01*dim(rfvotes)[1])))
      
      for(i in 1:length(e$Var1)){
        badinputs <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% select(inputnames) %>% c()
        badinputs.m <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% summarize(mean=mean(bestvote)) %>% c()
        badinput_ids <- unique(c(badinput_ids, badinputs$inputnames))
        classBestscore <- rfvotes %>% as.tibble() %>%  filter(., label == e$Var1[i]) %>% summarize(median=median(bestvote)) %>% c()
        #print(paste("The number of input dropped for",e$Var1[i],"class is", length(badinputs$inputnames),sep = " "))
        class_badinput_stats <- data.frame(class=e$Var1[i], classInputSize=e$Freq[i], allBestscoreMedian=Bestscore$median, classBestscoreMedian=classBestscore$median, tossedInput=length(badinputs$inputnames), tossedBestvoteMean=badinputs.m$mean, iteration=round_n)
        badinput_stats <- rbind(badinput_stats, class_badinput_stats)
      }
      
      badinput_stats[is.nan(badinput_stats)] <- 0
      toss_n <- badinput_stats %>% as.tibble() %>% filter(., iteration == round_n) %>% summarise(n=sum(tossedInput)) %>% c()
      toss_n <- toss_n$n
      
      print(badinput_stats)
      
      #filter input using the bad input list generated in the previous iteration
      trainingData <- trainingData[which(!rownames(trainingData) %in% badinput_ids ),]
      
      #run the RF again with the updated training set
      rf <- randomForest(CellType~., data = trainingData, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData$CellType)))
      print(rf)
      #Check again to see if there is room to improve:
      rfvotes <- as.data.frame(rf$votes)
      rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
      rfvotes$label <- apply(rf$votes, 1, function(x)  names(x)[which.max(x)] )
      rfvotes$inputnames <- rownames(rfvotes)
      e <- as.data.frame(table(rfvotes$label))
      Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
      
      #update the round
      round_n = round_n + 1
    }#closes the while loop
    
  }#closes the improve option
  
  save(rf, file=paste(run.name,".RF_model.Robj",sep = ""))
  return(rf)
}

PlotPredictions <- function(SeuratObject, model, save.pdf=T, outputFilename="plotpredictions"){
  #Evaluate model prediction accuracy:
  conf.mat <- model$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>%
    mutate(freq = 100*value/sum(value))
  
  class_n = length(model$classes)
  
  pdf(paste(outputFilename,".pdf",sep=""),width= 1.5*class_n, height = 1.5*class_n)
  
  
  require(gridExtra)
  
  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types")+
    scale_x_discrete(name ="Cell Types Classes")
  
  errorSize <- as.data.frame(cbind(model$confusion[,"class.error"],
                                   head(colSums(model$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- 100-100*colSums(model$confusion)["class.error"]/length(head(colSums(model$confusion),-1))
  names(acc) <- "accuracy"
  
  p2 <- ggplot(errorSize)  + 
    geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3")+
    geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
    scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize),name="Class % Error rate (Dots)"))+
    labs(y="Class Size in Training (Bars)",title=paste("Model Prediction Accuracy is ",round(acc,digits = 2),"%", sep=""))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+
    scale_x_discrete(name ="Cell Type Classes")
  
  grid.arrange(p1, p2, nrow=2)
  
  #Prediction outputs
  FeaturePlot(object = SeuratObject, 
              features.plot = model$classes, 
              cols.use = c("grey", "blue"), 
              reduction.use = "tsne")
  
  TSNEPlot(SeuratObject, group.by="Prediction",do.label=T)
  
  FeaturePlot(SeuratObject, features.plot = "BestVotesPercent", no.legend = F, cols.use = c("gray","red"))
  
  FeaturePlot(SeuratObject, features.plot = "KLe", no.legend = F, cols.use = c("gray","purple"))
  
  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=Prediction,fill=Prediction),stat = "count")+
    geom_violin(aes(x=Prediction,y=BestVotesPercent*max(table(SeuratObject@meta.data$Prediction)),fill=Prediction))+ 
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$Prediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))
  
  p4 <- ggplot(SeuratObject@meta.data, aes(KLe, Diff, color= PredictionStatus))+ geom_point(size=.6)
  
  grid.arrange(p3, p4, nrow=2)
  
  dev.off()
}

CellTyper <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions"){
  
  library(caret)
  library(randomForest)
  library(tidyverse)
  
  if(!missing(SeuratObject)){
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)
  colnames(testExpSet) <- make.names(colnames(testExpSet))
  #Prepare Test Expression set
  testsub <- testExpSet[,which(colnames(testExpSet) %in% attributes(model$terms)$term.labels)]
  missingGenes <- attributes(model$terms)$term.labels[which(!attributes(model$terms)$term.labels %in% colnames(testExpSet))]
  print(missingGenes)
  missingGenes.df <- data.frame(matrix(0, ncol = length(missingGenes), nrow = length(rownames(testExpSet))))
  colnames(missingGenes.df) <- missingGenes
  TestData <- cbind(testsub, missingGenes.df)
  TestData <- TestData[,attributes(model$terms)$term.labels]
  cat("Number of Features (genes) to be considered is", length(colnames(testsub)), '\n', "Number of missing Features set to zero is", length(missingGenes), '\n', sep = ' ')
  
  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class_n <- length(model$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class_n, class_n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="response")
  
  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class_n)  
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class_n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class_n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class_n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()
  
  if(missing(priorLabels)){
    print("Prior class labels are not provided!")
    
  }else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    testPred <- cbind(testPred, priorLabels)
    
    #Plot the crosscheck here:
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    crx <- testPred %>% group_by(Prior, Intermediate, Prediction) %>% tally() %>% as.data.frame()
    
    p5 <- ggplot(crx,aes(y = n, axis1 = Prior, axis2 = Intermediate, axis3 = Prediction )) +
      geom_alluvium(aes(fill = Prediction), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Prior", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")
    
    cowplot::save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class_n, base_width = 1.2*class_n)
    
  }
  
  if(!missing(SeuratObject)){
    
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)
    
    PlotPredictions(SeuratObject = SeuratObject, model = model, outputFilename = outputFilename)
    
    return(SeuratObject)
  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function

