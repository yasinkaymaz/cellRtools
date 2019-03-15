#Test functions


DFT <- function(tree,node,path = NULL,include.children = FALSE,only.children = FALSE) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if(! only.children){
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <-c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (! only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}
GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}
BuildRFClassifier <- function(object,training.genes = NULL,training.classes = NULL,verbose = TRUE,node=node,...) {
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
  classifier <- caret::train(CellType~., data=training.data, trControl=train_control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=50)
  return(classifier)
}
BuildClusterTree <- function(object, genes.use = NULL,  pcs.use = NULL, SNN.use = NULL, do.plot = TRUE, do.reorder = FALSE, reorder.numeric = FALSE, show.progress = TRUE) {
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  ident.names <- as.character(x = unique(x = object@ident))
  if (! is.null(x = genes.use)) {
    genes.use <- intersect(x = genes.use, y = rownames(x = object@data))
    data.avg <- AverageExpression(
      object = object,
      genes.use = genes.use,
      show.progress = show.progress
    )
    data.dist <- dist(t(x = data.avg[genes.use, ]))
  }
  if (! is.null(x = pcs.use)) {
    data.pca <- AveragePCA(object = object)
    data.dist <- dist(t(x = data.pca[pcs.use,]))
  }
  if (! is.null(x = SNN.use)) {
    num.clusters <- length(x = ident.names)
    data.dist <- matrix(data = 0, nrow = num.clusters, ncol = num.clusters)
    rownames(data.dist) <- ident.names
    colnames(data.dist) <- ident.names
    for (i in 1:(num.clusters - 1)) {
      for (j in (i + 1):num.clusters) {
        subSNN <- SNN.use[
          match(
            x = WhichCells(object = object, ident = ident.names[i]),
            table = colnames(x = SNN.use)
          ), # Row
          match(
            x = WhichCells(object = object, ident = ident.names[j]),
            table = rownames(x = SNN.use)
          ) # Column
          ]
        d <- mean(subSNN)
        if (is.na(x = d)) {
          data.dist[i, j] <- 0
        } else {
          data.dist[i, j] <- d
        }
      }
    }
    diag(x = data.dist) <- 1
    data.dist <- dist(data.dist)
  }
  data.tree <- as.phylo(x = hclust(d = data.dist))
  object@cluster.tree[[1]] <- data.tree
  if (do.reorder) {
    old.ident.order <- sort(x = unique(x = object@ident))
    data.tree <- object@cluster.tree[[1]]
    all.desc <- GetDescendants(tree = data.tree, node = (data.tree$Nnode + 2))
    all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 1)]]
    object@ident <- factor(x = object@ident, levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
      object <- SetIdent(
        object = object,
        cells.use = object@cell.names,
        ident.use = as.integer(x = object@ident)
      )
      object@meta.data[object@cell.names, "tree.ident"] <- as.integer(x = object@ident)
    }
    object <- BuildClusterTree(
      object = object,
      genes.use = genes.use,
      pcs.use = pcs.use,
      do.plot = FALSE,
      do.reorder = FALSE,
      show.progress = show.progress
    )
  }
  if (do.plot) {
    PlotClusterTree(object)
  }
  return(object)
}

PlotClusterTree <- function(object, ...) {
  if (length(x = object@cluster.tree) == 0) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- object@cluster.tree[[1]]
  plot.phylo(x = data.tree, direction = "downwards", ...)
  nodelabels()
}

PlotCellPredictionTree <- function(tree){

  ape::plot.phylo(tree,direction = "downwards")
  dd <- data.frame("L"=c(.8,.9,.5),"R"=c(.2,.1,.5))
  ape::nodelabels(thermo = dd, horiz = TRUE,cex=1.2)
}
  


SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}
PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- pkgs %in% rownames(x = installed.packages())
  if (error && any(!package.installed)) {
    stop(
      "Cannot find ",
      paste(pkgs[!package.installed], collapse = ', '),
      "; please install"
    )
  }
  invisible(x = package.installed)
}
AssessSplit <- function(object,node,cluster1,cluster2,genes.training = NULL,print.output = TRUE,...) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@cluster.tree[[1]]
  if (!missing(x = node)) {
    if (!missing(x = cluster1) || !missing(x = cluster2)) {
      warning("Both node and cluster IDs provided. Defaulting to using node ID")
    }
    possible.nodes <- c(
      DFT(tree = tree, node = tree$edge[,1][1]),
      tree$edge[,1][1]
    )
    if (!node %in% possible.nodes) {
      stop("Not a valid node")
    }
    split <- tree$edge[which(x = tree$edge[,1] == node), ][,2]
    group1 <- DFT(tree = tree, node = split[1], only.children = TRUE)
    group2 <- DFT(tree = tree, node = split[2], only.children = TRUE)
    if (any(is.na(x = group1))) {
      group1 <- split[1]
    }
    if (any(is.na(x = group2))) {
      group2 <- split[2]
    }
  } else {
    group1 <- cluster1
    group2 <- cluster2
  }
  group1.cells <- WhichCells(object = object, ident = group1)
  group2.cells <- WhichCells(object = object, ident = group2)
  assess.data <- SubsetData(
    object = object,
    cells.use = c(group1.cells, group2.cells)
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells.use = group1.cells,
    ident.use = "L"
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells.use = group2.cells,
    ident.use = "R"
  )
  rfc <- BuildRFClassifier(
    object = assess.data,
    training.genes = genes.training,
    training.classes = assess.data@ident,
    node=node,
    ...
  )

  oobe <- rfc$finalModel$confusion[,"class.error"]
  if (print.output) {
    message(paste0("Out of Bag Error::: ", round(x = oobe, digits = 4) * 100, "%"))
  }

  return(rfc)
}
AssessNodes <- function(object,node.list,all.below = FALSE,genes.training = NULL) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@cluster.tree[[1]]
  if (missing(x = node.list)) {
    node.list <- GetAllInternalNodes(tree = tree)
  } else {
    possible.nodes <- GetAllInternalNodes(tree = tree)
    if (any(!node.list %in% possible.nodes)) {
      stop(paste(
        node.list[!(node.list %in% possible.nodes)],
        "not valid internal nodes"
      ))
    }
    if (length(x = node.list == 1) && all.below) {
      node.list <- c(node.list, DFT(tree = tree, node = node.list))
    }
  }
  rfs <- vector("list", length = max(node.list))
  for(i in node.list){
    rfs[[i]] <- AssessSplit(object = object, node=i, genes.training = genes.training, print.output = F,verbose=F)
  }
  names(rfs) <- seq_along(rfs)
  rfs[sapply(rfs, is.null)] <- NULL

  return( rfs)
}
GetAncestorsPath <- function(tree, leafNode){
  ancestors <- c(); leafside <- c();
  ni = leafNode;
  parent <- tree$edge[which(x = tree$edge[,2] == leafNode), ][1]
  while(!is.na(parent)){
    if(tree$edge[which(x = tree$edge[,1] == parent), ][,2][1] == ni){
      leafside <- c(leafside, "L")
    }else{
      leafside <- c(leafside, "R")
    }
    ancestors <- c(ancestors, parent)
    ni <- parent
    parent <- tree$edge[which(x = tree$edge[,2] == parent), ][1]
  }
  path <- paste(ancestors,leafside,sep = "")
  map <- list(leafNode, path)
  return(map)
}
ClassifierTreeTraverser <- function(testExpSet,tree, CLoc.list){

  Htable <- data.frame(cells = rownames(testExpSet))
  Dtable <- data.frame(cells = rownames(testExpSet))
  internalNodes <- GetAllInternalNodes(tree = tree)

  for (ni in internalNodes){

    modelname <- ni
    model <- CLoc.list[[as.character(ni)]]
    print(paste("Predicting with local classifier (CLoc.i) at node ",modelname,"...",sep = " "))
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

    #if((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )){
    #  warning("A significant portion of features are missing...")
    #}

    rm(testsub, missingGenes, missingGenes.df)
    gc()
    #Predict
    library(entropy)
    testPred <- as.data.frame(predict(model, TestData, type = "prob"))
    colnames(testPred) <- paste(ni, colnames(testPred), sep = "")

    chosenGate <- as.data.frame(predict(model, TestData, type = "raw"))
    colnames(chosenGate) <- as.character(ni)

    Htable <- cbind(Htable, testPred)
    Dtable <- cbind(Dtable, chosenGate)
  }#closes models for loop

  print(head(Dtable))
  print(class(Dtable))

  checknodes <- function(x, internalNodes, tree) {
    ni = internalNodes[1];
    while (ni %in% internalNodes){
      if (x[as.character(ni)] == "L"){
        ni = tree$edge[which(x=tree$edge[,1] == ni  ),2][1]
      }else{
        ni = tree$edge[which(x=tree$edge[,1] == ni  ),2][2]
      }
    }
    nt = ni #terminal node
    return(nt)
  }

  Dtable$Decision <- apply(Dtable, 1, function(x) checknodes(x, internalNodes=internalNodes, tree=tree))
  print(Dtable)
  CTTtables <- list(Htable, Dtable)
  return(CTTtables)
}
ClassProbCalculator <- function(tree, Htable){
  #CTip_table <- data.frame(matrix(ncol = 0, nrow = length(rownames(Htable))))
  CTip_table <- data.frame(cells = rownames(Htable))
  for(nt in tree$tip.label){
    map <- GetAncestorsPath(tree = tree, leafNode = nt)

    nt_prob <- data.frame(matrixStats::rowProds(as.matrix(Htable[,map[[2]]])))
    colnames(nt_prob) <- paste(nt, "classProb",sep = "_")
    print(head(nt_prob))
    CTip_table <- cbind(CTip_table, nt_prob)
  }
  print(head(CTip_table))
  return(CTip_table)
}
HTyper2 <- function(SeuratObject, tree, testExpSet, models, priorLabels, outputFilename="plotpredictions"){

  #models is a list of of rf models
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

  CTTtables <- ClassifierTreeTraverser(testExpSet = testExpSet, tree = tree, CLoc.list = CLoc.list)

  Ctable <- ClassProbCalculator(tree = tree, Htable = CTTtables[[1]] )
  Ctable$HRFPrediction <- str_remove(colnames(Ctable)[apply(Ctable,1,which.max)],"_classProb")

  if(!missing(SeuratObject)){
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(cbind(Ctable,CTTtables[[2]])))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, Ctable, CTTtables[[2]])

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function
CellTyperTrainer <- function(ExpressionData, CellLabels, run.name, do.splitTest=F, PCs, improve=T){
  library(randomForest)
  #library(rfUtilities)
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
    trainingData <- prepareDataset(ExpressionData = ExpressionData, CellLabels = CellLabels, PCs = PCs, run.name = run.name)
  }

  #Added: "sampsize=c(table(trainingData$CellType))". Revisit this later to make sure it is working as expected...
  rf <- randomForest(CellType~., data = trainingData, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData$CellType)))
  print(rf)
  save(rf, file=paste(run.name,".RF_model_notImproved.Robj", sep = ""))

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
    #while (Bestscore$median >= currentscore && toss_n > round(0.01*dim(rfvotes)[1]) ){
    while (round_n < 2 ){#run this only once...

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
    crx <- testPred %>% group_by(Prior, Intermediate, Prediction) %>% tally() %>% as.data.frame()

    p5 <- ggplot(crx,aes(y = n, axis1 = Prior, axis2 = Intermediate, axis3 = Prediction )) +
      geom_alluvium(aes(fill = Prediction), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Prior", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    cowplot::save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 16, base_width = 20)

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
