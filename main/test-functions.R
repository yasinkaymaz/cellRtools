#Test functions

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
    ni <- internalNodes[1]
    while (ni %in% internalNodes){
      if (x[as.character(ni)] == "L"){
        ni = tree$edge[which(x=tree$edge[,1] == ni  ),2][1]
      }else{
        ni = tree$edge[which(x=tree$edge[,1] == ni  ),2][2]
      }
    }
    nt <- ni #terminal node
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
    filter_p <- 0.05
    bestvote_cutoff <- 0.7
    badinput_ids <- NULL
    round_n <- 1
    badinput_stats <- data.frame()
    currentscore <- Bestscore$median
    toss_n <- dim(rfvotes)[1]

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
      round_n <- round_n + 1
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

CrossCheck <- function(PriorPostTable,outputprefix=""){
  #This function takes a table with two columns of which first is for prior cell labels and second is for predicted cell class.
  library(tidyverse)
  library(alluvial)
  library(ggalluvial)

  crx <- PriorPostTable %>% group_by_at(vars(one_of(names(PriorPostTable))))%>% tally() %>% arrange(desc(n))

  p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2] )) +
    geom_alluvium(aes_string(fill = names(crx)[2]), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "red") +
    geom_label(stat = "stratum", label.strata = TRUE) +
    scale_x_discrete(limits = c("PriorLabels","FinalPrediction"), expand = c(.05, .05)) +
    ggtitle("Predictions Cross-Check")
  pdf(paste(outputprefix,".prediction-crosscheck.pdf",sep=""),width = 20,height = 15)
  print(p5)
  dev.off()

}

TipBias <- function(tree, confmat){
  #confmat is table() of Prior/Prediction comparison
  err.rate <- NULL
  for(i in 1:length(tree$tip.label)){
    nn <- length(GetAncestorsPath(tree=tree, i)[[2]])
    leafErr <- 1-confmat[tree$tip.label[i],tree$tip.label[i]]/sum(confmat[tree$tip.label[i],])
    print(paste(i,nn,leafErr,sep = "    "))
    err.rate <- c(err.rate, leafErr)
  }
  err.rate <- data.frame(err.rate)
  rownames(err.rate) <- tree$tip.label
  return(err.rate)
}

#' A hierarchical Classifier Tree Traverser function.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param hiemods models list from HieRandForest function.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
CTTraverser <- function(Query, tree, hiemods, thread=NULL){

  node.list <- DigestTree(tree = tree)
  #Create a table for storing node probabilities.
  Scores <- data.frame(row.names = colnames(Query))
  Pvotes <- data.frame(row.names = colnames(Query))
  QueWs <- data.frame(row.names = colnames(Query))
  QueCers <- data.frame(row.names = colnames(Query))
  #fi=node.list[length(node.list)]+1
  #node.list <- c(node.list, fi)
  #Query_rand <- RandomizeR(df = t(Query), n = 5)

  if(is.null(thread)){
    for(i in node.list){
      #nodeModel <- hiemods[[as.character(i)]][[1]]
      nodeModel <- hiemods[[as.character(i)]]
      #c_f <- length(nodeModel$levels)
      #Create QueData:
      P_dicts <- FixLab(xstring = nodeModel$finalModel$xNames)
      nodeQueData <- Query[which(rownames(Query) %in% P_dicts), ]
      nodeQueData <- t(nodeQueData)
      #Add the missing features matrix
      mp <- P_dicts[which(!P_dicts %in% colnames(nodeQueData))]
      mp_df <- data.frame(matrix(0,
                                 ncol = length(mp),
                                 nrow = length(colnames(Query))))
      colnames(mp_df) <- mp
      nodeQueData <- cbind(nodeQueData, mp_df)

      #Create randomized node Query data:
      #nodeQueData_rand <- Query_rand[, which(colnames(Query_rand) %in% P_dicts)]
      #nodeQueData_rand <- cbind(nodeQueData_rand, mp_df)
      #nodeQueWs <- graWeighteR(model = nodeModel, QueData = nodeQueData_rand )

      #Calculate the probability weights of each class by random permutation:
      nodeQueWs <- graWeighteR(model = nodeModel, QueData = nodeQueData )
      #nodeQueWs2 <- graWeighteR2(model = nodeModel, QueData = nodeQueData )

      #Tally votes for class from the local model:
      nodePvotes <- PvoteR(model = nodeModel, QueData = nodeQueData)

      #Estimate Certainty of prediction probabilities per class:
      nodeQueCers <- ceR(qP = nodePvotes, qW = nodeQueWs)
      #Calculate node Scores:
      #nodeScores <- scoR(model = nodeModel,
      #                   format = "prob",
      #                   QueData = nodeQueData,
      #                   node = i)
      nodeScores <- nodePvotes/nodeQueWs
      #nodeScores <- nodePvotes
      colnames(nodeScores) <- paste(i, colnames(nodeScores), sep = "")
      Scores <- cbind(Scores, nodeScores)
      #if(i != fi){
        Pvotes <- cbind(Pvotes, nodePvotes)
        QueWs <- cbind(QueWs, nodeQueWs)
        QueCers <- cbind(QueCers, nodeQueCers)
      #}
    } #closes the for loop.
    S_path_prod <- ClassProbCalculator(tree = tree, nodes_P_all = Scores)
    #S_path_prod <- ClassProbCalculator2(tree = tree, nodes_P_all = Scores)
    HieMetrxObj <- new(Class = "HieMetrics",
                       Pvotes = Pvotes,
                       QueWs = QueWs,
                       QueCers = QueCers,
                       Scores = S_path_prod)

  }else{# FIX THIS PART! thread is specified. For now, use this only when running on bigMem machines.
    library(doParallel)
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    nodeProb <- foreach(i=node.list, .inorder = TRUE, .combine=cbind) %dopar% {
      scoR(model = hiemods[[as.character(i)]],
           format = "prob",
           QueData = Query,
           node = i)
    }
    stopCluster(cl)
    ProbTab <- cbind(ProbTab, nodeProb)
  }

  return(HieMetrxObj)
}


#' A hierarchical Classifier Tree Traverser function.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param hiemods models list from HieRandForest function.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
CTTraverser4 <- function(Query, tree, hiemods, thread=NULL){

  node.list <- DigestTree(tree = tree)
  #Create a table for storing node probabilities.
  Scores <- data.frame(row.names = colnames(Query))
  Pvotes <- data.frame(row.names = colnames(Query))
  QueWs <- data.frame(row.names = colnames(Query))
  QueCers <- data.frame(row.names = colnames(Query))

  if(is.null(thread)){
    for(i in node.list){
      nodeModel <- hiemods[[as.character(i)]]
      #c_f <- length(nodeModel$levels)
      #Create QueData:
      #P_dicts <- FixLab(xstring = nodeModel$finalModel$xNames)
      P_dicts <- FixLab(xstring = rownames(nodeModel$trainPCA$rotation))
      nodeQueData <- Query[which(rownames(Query) %in% P_dicts), ]
      nodeQueData <- t(nodeQueData)
      #Add the missing features matrix
      mp <- P_dicts[which(!P_dicts %in% colnames(nodeQueData))]
      mp_df <- data.frame(matrix(0,
                                 ncol = length(mp),
                                 nrow = length(colnames(Query))))
      colnames(mp_df) <- mp
      nodeQueData <- cbind(nodeQueData, mp_df)

      #Calculate the probability weights of each class by random permutation:
      nodeQueWs <- graWeighteR4(model = nodeModel, QueData = nodeQueData )
      #Tally votes for class from the local model:
      nodePvotes <- PvoteR4(model = nodeModel, QueData = nodeQueData)
      #Estimate Certainty of prediction probabilities per class:
      nodeQueCers <- ceR(qP = nodePvotes, qW = nodeQueWs)
      #Calculate node Scores:
      nodeScores <- nodePvotes/nodeQueWs
      colnames(nodeScores) <- paste(i, colnames(nodeScores), sep = "")
      Scores <- cbind(Scores, nodeScores)

      Pvotes <- cbind(Pvotes, nodePvotes)
      QueWs <- cbind(QueWs, nodeQueWs)
      QueCers <- cbind(QueCers, nodeQueCers)

    } #closes the for loop.
    S_path_prod <- ClassProbCalculator(tree = tree, nodes_P_all = Scores)
    HieMetrxObj <- new(Class = "HieMetrics",
                       Pvotes = Pvotes,
                       QueWs = QueWs,
                       QueCers = QueCers,
                       Scores = S_path_prod)

  }else{
    library(doParallel)
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    nodeProb <- foreach(i=node.list, .inorder = TRUE, .combine=cbind) %dopar% {
      scoR(model = hiemods[[as.character(i)]],
           format = "prob",
           QueData = Query,
           node = i)
    }
    stopCluster(cl)
    ProbTab <- cbind(ProbTab, nodeProb)
  }

  return(HieMetrxObj)
}


#' A function to tally votes from the model locally trained.
#' @param model
#' @param QueData is the prepared data matrix ready to be used in predict
#' @param format type of prediction output, "prob" or "resp".
PvoteR4 <- function(model, QueData, format="prob", node=NULL){

  QueDataPCs <- predict(model$trainPCA, newdata = QueData)
  QueDataPCs <- as.data.frame(QueDataPCs)

  # if(is.null(node)){
  #   QuePvotes <- as.data.frame(predict(model, newdata=QueDataPCs, type = "prob", scale=T, center=T))
  # } else{
  #   if(format == "prob"){
  #     QuePvotes <- as.data.frame(predict(model, newdata=QueDataPCs, type = "prob", scale=T, center=T))
  #     colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
  #   } else {
  #     QuePvotes <- as.data.frame(predict(model, newdata=QueDataPCs, type = "raw", scale=T, center=T))
  #     colnames(QuePvotes) <- as.character(node)
  #   }
  # }
  if(is.null(node)){
    QuePvotes <- as.data.frame(predict(model, newdata=QueDataPCs, type = "prob", scale=F, center=F))
  } else{
    if(format == "prob"){
      QuePvotes <- as.data.frame(predict(model, newdata=QueDataPCs, type = "prob", scale=F, center=F))
      colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
    } else {
      QuePvotes <- as.data.frame(predict(model, newdata=QueDataPCs, type = "raw", scale=F, center=F))
      colnames(QuePvotes) <- as.character(node)
    }
  }
  return(QuePvotes)
}

#' A function to calculate gravitity center (weight center) of the probability distributions
#' among classes of the node. By permutation.
#' @param model of the node.
#' @param QueData same matrix used in PvoteR().
#' @return QueWeights a set of probability weightes per class to be used in asymetrix entropy estimations.
graWeighteR4 <- function(model, QueData){
  #Randomizing only feature space
  #QueData_R <- RandomizeR(df = QueData, n = 20)
  QueData_R <- Shuffler(df = QueData)
  #QueData_R <- Rander(df = QueData)
  pvts_R <- PvoteR4(model = model, QueData = QueData_R)
  #Ws <- apply(pvts_R, 2, mean) + apply(pvts_R, 2, sd)
  #Ws <- apply(pvts_R, 2, mean) + apply(pvts_R, 2, sd)/sqrt(dim(pvts_R)[1])
  Ws <- apply(pvts_R, 2, mean)
  #Ws <- colMeans(PvoteR(model = model, QueData = QueData_R))
  QueWeights <- t(as.data.frame(Ws))[rep(1, each=nrow(QueData)), ]
  QueWeights <- as.data.frame(QueWeights)
  return(QueWeights)
}


#' A function to tally votes from the model locally trained.
#' @param model
#' @param QueData is the prepared data matrix ready to be used in predict
#' @param format type of prediction output, "prob" or "resp".
CalibratedPvoteR.LR <- function(model, refmod, QueData, format="prob", node=NULL){
  if(is.null(node)){
    QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
  } else{
    if(format == "prob"){
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
      colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
    } else {
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "raw", scale=T, center=T))
      colnames(QuePvotes) <- as.character(node)
    }
  }
  #Calibrate the probabilities:
  QuePvotes.cal <- data.frame(row.names = rownames(QuePvotes))
  for(class in model$finalModel$classes){
    prep.cal <- predict(refmod@sigmoids[[class]],
                        data.frame(x=QuePvotes[, class]),
                        type="response")
    Calibrated_and_Scaled = (prep.cal - min(prep.cal))/(max(prep.cal) - min(prep.cal))
    QuePvotes.cal <- cbind(QuePvotes.cal, data.frame(class=Calibrated_and_Scaled))
  }
  colnames(QuePvotes.cal) <- model$finalModel$classes

  return(QuePvotes.cal)
}

#' A function to calculate gravitity center (weight center) of the probability distributions
#' among classes of the node. By permutation.
#' @param model of the node.
#' @param QueData same matrix used in PvoteR().
#' @return QueWeights a set of probability weightes per class to be used in asymetrix entropy estimations.
graWeighteR <- function(model, QueData){
  #Randomizing only feature space
  #QueData_R <- RandomizeR(df = QueData, n = 20)
  #QueData_R <- Rander(df = QueData)
  QueData_R <- Shuffler(df = QueData)
  pvts_R <- PvoteR(model = model, QueData = QueData_R)
  #Ws <- apply(pvts_R, 2, mean) + apply(pvts_R, 2, sd)
  #Ws <- apply(pvts_R, 2, mean) + apply(pvts_R, 2, sd)/sqrt(dim(pvts_R)[1])
  Ws <- apply(pvts_R, 2, mean)
  #Ws <- colMeans(PvoteR(model = model, QueData = QueData_R))
  QueWeights <- t(as.data.frame(Ws))[rep(1, each=nrow(QueData)), ]
  QueWeights <- as.data.frame(QueWeights)
  return(QueWeights)
}

graWeighteR3 <- function(model, QueData){
  allZero <- QueData
  allZero[1:ncol(allZero),] <- 0
  allZero[,1:ncol(allZero)] <- 0
  pvts_R <- PvoteR(model = model, QueData = allZero)
  Ws <- apply(pvts_R, 2, mean)
  #Ws <- colMeans(PvoteR(model = model, QueData = QueData_R))
  QueWeights <- t(as.data.frame(Ws))[rep(1, each=nrow(QueData)), ]
  QueWeights <- as.data.frame(QueWeights)
  return(QueWeights)
  return(QueWeights)
}

graWeighteR2 <- function(model, QueData){
  #Randomizing only feature space
  QueData_R <- RandomizeR(df = QueData, n = 20)
  QueData_R <- QueData_R[sample(rownames(QueData_R), size = nrow(QueData)),]
  pvts_R <- PvoteR(model = model, QueData = QueData_R)
  pvts_R <- as.data.frame(pvts_R)
  return(pvts_R)
}



#' The predictor function. REDUNDANT!!!
#' @param model
#' @param QueData
#' @param format type of prediction output, "prob" or "resp".
scoR <- function(model, QueData, format="prob", node=NULL){

  if(is.null(node)){
    QuePred <- as.data.frame(predict(model, QueData, type = "prob", scale=T, center=T))
  } else{
    if(format == "prob"){
      c_f <- length(model$levels) #correction factor; class size
      QuePred <- as.data.frame(predict(model, QueData, type = "prob", scale=T, center=T))
      #      QuePred <- as.data.frame(t(apply(QuePred, 1, function(x) KLeCalc(x)*x)))
      QuePred <- QuePred*c_f
      colnames(QuePred) <- paste(node, colnames(QuePred), sep = "")
    } else {
      QuePred <- as.data.frame(predict(model, QueData, type = "raw", scale=T, center=T))
      colnames(QuePred) <- as.character(node)
    }
  }
  return(QuePred)
}



ClassProbCalculator2 <- function(tree, nodes_P_all){
  P_path_prod <- data.frame(row.names = rownames(nodes_P_all))
  clabs <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  fi=length(c(tree$tip.label, tree$node.label))+1
  for(cl in clabs){
    map <- GetAncestPath(tree = tree, class = cl)
    if(cl %in% tree$tip.label){
      map <- c(map, paste(fi, cl, sep = ""))
    }
    if(length(map) > 1){
      nt_prob <- data.frame(apply(nodes_P_all[, map], 1, prod))
    }else{
      nt_prob <- data.frame(nodes_P_all[, map])
    }
    colnames(nt_prob) <- cl
    P_path_prod <- cbind(P_path_prod, nt_prob)
  }
  return(P_path_prod)
}



#' A function for evalating the uncertainty.
#' @param ScoreObs P_path_prod for observed scores
#' @param ProbCert Certainty scores.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
ScoreEval <- function(ScoreObs, ProbCert, tree, alpha=0){

  df <- data.frame(row.names = rownames(ScoreObs))
  for(i in 1:length(ProbCert[,1])){
    #candits <- colnames(ProbCert)[ProbCert[i,] > alpha]
    #candits <- CandidateDetector(PCertVector = ProbCert.logic[i,], tree=tree)
    candits <- CandidateDetector2(PCertVector = ProbCert[i,], tree=tree, alpha=alpha)
    #candits <- CandidateDetector3(PCertVector = ProbCert[i,], tree=tree, alphas=alpha)
    if(length(candits) == 0){
      classL <- "Undetermined"
      classU <- max(ProbCert[i,])
      classS <- max(ScoreObs[i,])
    }else if(length(candits) == 1){
      classL <- candits
      #classU <- ProbCert[i,classL]
      AncPath <- GetAncestPath(tree = tree, class = classL, labels = T)
      if(length(AncPath) > 1){classU <- rowMeans(ProbCert[i, AncPath])}else{classU <- ProbCert[i, AncPath]}
      classS <- ScoreObs[i, classL]
    }else{
      classL <- colnames(ScoreObs[i,candits])[apply(ScoreObs[i,candits], 1, which.max)]
      #classL <- colnames(ProbCert[i,candits])[apply(ProbCert[i,candits], 1, which.max)]
      #classU <- ProbCert[i,classL]
      AncPath <- GetAncestPath(tree = tree, class = classL, labels = T)
      if(length(AncPath) > 1){classU <- rowMeans(ProbCert[i, AncPath])}else{classU <- ProbCert[i, AncPath]}
      classS <- ScoreObs[i, classL]
    }

    df <- rbind(df, data.frame(row.names = rownames(ProbCert[i,]),
                               Score = classS,
                               Certainty = classU,
                               Projection = classL))
  }
  return(df)
}



#' A function to find the full path of all true ancestors:
#' @param PCertVector Probability certainty table.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
CandidateDetector <- function(PCertVector, tree){

  labs_l <- c(tree$tip.label, tree$node.label)#The order is important! tips first. Don't change!#New
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"] #Double-check#New
  Path_nodes_of_candits <- NULL
  CandidNodes <- NULL
  for(node.lab in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    if(all(as.logical(PCertVector[AncPath])) & !(AncPath[1] %in% Path_nodes_of_candits)){
      CandidNodes <- c(CandidNodes, AncPath[1])
      Path_nodes_of_candits <- unique(c(Path_nodes_of_candits, AncPath[2:length(AncPath)]))
    }
  }
  return(CandidNodes[!CandidNodes %in% Path_nodes_of_candits])
}

#' A function to find the full path of all true ancestors:
#' @param PCertVector Probability certainty table.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
CandidateDetector2 <- function(PCertVector, tree, alpha=0){
  labs_l <- c(tree$tip.label, tree$node.label)#The order is important! tips first. Don't change!#New
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"] #Double-check#New
  Path_nodes_of_candits <- NULL
  CandidNodes <- NULL
  for(node.lab in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    if( (mean(as.numeric(PCertVector[AncPath])) > alpha) &
        !(AncPath[1] %in% Path_nodes_of_candits) #&
        #as.numeric(PCertVector[AncPath[1]]) > alpha
        ){
      CandidNodes <- c(CandidNodes, AncPath[1])
      Path_nodes_of_candits <- unique(c(Path_nodes_of_candits, AncPath[2:length(AncPath)]))
    }
  }
  return(CandidNodes)
}

#' A function to find the full path of all true ancestors:
#' @param PCertVector Probability certainty table.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param alphas a list of numeric values of node specific alpha thresholds. Learned from reference data.
CandidateDetector3 <- function(PCertVector, tree, alphas){
  labs_l <- c(tree$tip.label, tree$node.label)#The order is important! tips first. Don't change!#New
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"] #Double-check#New
  Path_nodes_of_candits <- NULL
  CandidNodes <- NULL
  for(node.lab in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    #if( all(PCertVector[AncPath] > alphas[[node.lab]]) &
    if( (mean(as.numeric(PCertVector[AncPath])) > mean(as.numeric(alphas[[node.lab]]))) &
        !(AncPath[1] %in% Path_nodes_of_candits)
    ){
      CandidNodes <- c(CandidNodes, AncPath[1])
      Path_nodes_of_candits <- unique(c(Path_nodes_of_candits, AncPath[2:length(AncPath)]))
    }
  }
  return(CandidNodes)
}


SubsetTData2 <- function(Tdata, tree, node){

  labs_l <- c(tree$tip.label, tree$node.label)

  if(labs_l[node] %in% tree$tip.label){
    childNodes <- list(labs_l[node])
    SubTdata <- droplevels(Tdata[which(Tdata$ClassLabels %in% childNodes[[1]]), ])

  }else{

  # 1. Extract the data under the node. Subsample if necessary.
  childNodes <- GetChildNodeLeafs(tree = tree, node = node)
  SubTdata <- NULL
  #loop through child nodes that are not null in the list.
  for (i in which(lapply(childNodes, length) > 0)){
    Subdata <- droplevels(Tdata[which(Tdata$ClassLabels %in% childNodes[i][[1]]), ])
    if (i > length(tree$tip.label)){# if the node is not a leaf node, then
      if(!is.null(tree$node.label)){# if labels for subnodes exist
        labels <- c(tree$tip.label, tree$node.label)
        #Replace labels with subnode labels.
        Subdata$ClassLabels <- as.factor(labels[i])
      } else {
        #if subnodes don't exist, replace class tip labels with childnode label.
        Subdata$ClassLabels <- as.factor(i)
      }
    } else {#if the node is a terminal leaf
      #Replace class tip labels with Leaf labels.
      Subdata$ClassLabels <- as.factor(childNodes[i][[1]])
    }
    #Combine with all other child node data
    SubTdata <- rbind(SubTdata, Subdata)
  }
 }
  return(SubTdata)
}




#' A function to create a local classifier for a given node.
#' @param Rdata
#' @param tree
#' @param node
#' @param f_n number of features to be included in local classifier.
NodeTrainer <- function(Rdata, tree, node, f_n=200, tree_n=500, ...){

  node.Data <- SubsetTData2(Tdata = Rdata, tree = tree, node = node)
  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  node.Data <- node.Data[, apply(node.Data, 2, var) != 0]
  #First select the highly variable genes that correlate with the PCs
  P_dict <- FeatureSelector(Data = node.Data,
                             ClassLabels = node.ClassLabels,
                             num = 2000,
                             ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  #Then, select the genes as predictors if statistically DE between the classes.
  P_dict <- FeatureSelector2(Data = node.Data,
                            ClassLabels = node.ClassLabels,
                            num = f_n,
                            ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.internal? is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)

  if(labs_l[node] != "TaxaRoot"){
  parent.t <- tree$edge[which(x = tree$edge[, 2] == node), ][1]
  for(n in 1:length(labs_l)){
    if(labs_l[n] != "TaxaRoot"){
      parent.p <- tree$edge[which(x = tree$edge[, 2] == n), ][1]
      if(parent.t == parent.p & node != n){
        print(paste("Training for", labs_l[node], "and its sibling is:", labs_l[n]))
        #SubsetTData
        node.outData <- SubsetTData2(Tdata = Rdata, tree = tree, node = n)
        node.outData <- droplevels(subset(node.outData, select=c(P_dict)))
        node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
        if(dim(node.outData)[1] > 500){
          node.outData <- node.outData[sample(rownames(node.outData), size = 500), ]
        }
        node.Data <- rbind(node.Data, node.outData)
       }
     }
   }
  }
  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)
  node.mod <- caret::train(ClassLabels~.,
                           data = node.Data,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = FALSE,
                           proximity = FALSE,
                           outscale = FALSE,
                           preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)

  return(node.mod)
}

NodeTrainer3 <- function(Rdata, tree, node, f_n=200, tree_n=500, switchBox='off', ...){

  node.Data <- SubsetTData(Tdata = Rdata, tree = tree, node = node)
  #node.ClassLabels <- node.Data$ClassLabels
  #node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  #node.Data <- node.Data[, apply(node.Data, 2, var) != 0]

  #node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)
  if(labs_l[node] != "TaxaRoot"){
    childNodes <- GetChildNodeLeafs(tree = tree, node = node)
    childLeafs <- NULL
    for (i in which(lapply(childNodes, length) > 0)){
      childLeafs <- c(childLeafs, childNodes[i][[1]])
    }
    outGroupLeafs <- tree$tip.label[!tree$tip.label %in% childLeafs]
    node.outData <- droplevels(Rdata[which(Rdata$ClassLabels %in% outGroupLeafs), ])
    if(dim(node.outData)[1] > 500){
      node.outData <- node.outData[sample(rownames(node.outData), size = 500), ]#500 can be replaced with a better #
    }
    node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
  }else{
    #Create OurGroup for the TaxaRoot
    node.outData <- droplevels(Rdata[sample(rownames(Rdata), size = 500), ])
    node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
    node.outData <- droplevels(subset(node.outData, select=-c(ClassLabels)))
    node.outData <- Shuffler(df = node.outData)
  }
  node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
  node.Data <- rbind(node.Data, node.outData)

  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))

  node.Data <- node.Data[, apply(node.Data, 2, var) != 0]
#$$$
  #Select top highly varible genes:
  vmr <- Seurat::FindVariableFeatures(t(expm1(node.Data)),verbose=F)
  vmr <- vmr[order(vmr$vst.variance.standardized, decreasing = T),]
  P_dict <- head(rownames(vmr), 4000)

  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
#$$$
    #First select the highly variable genes that correlate with the PCs
    P_dict <- FeatureSelector(Data = node.Data,
                              ClassLabels = node.ClassLabels,
                              num = 2000,
                              ...)
    node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
    #Then, select the genes as predictors if statistically DE between the classes.
    P_dict <- FeatureSelector2(Data = node.Data,
                               ClassLabels = node.ClassLabels,
                               num = f_n,
                               ...)

    #pool <- PairsPool(P_dict = P_dict, node.Data = node.Data, node.ClassLabels = node.ClassLabels)

    node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
    #node.Data <- droplevels(subset(node.Data, select=c(pool)))

  #node.Data$ClassLabels <- node.ClassLabels
  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)
  # node.mod <- caret::train(ClassLabels~.,
  #                          data = node.Data,
  #                          method = "rf",
  #                          norm.votes = TRUE,
  #                          importance = TRUE,
  #                          proximity = FALSE,
  #                          outscale = FALSE,
  #                          preProcess = c("center", "scale"),
  #                          ntree = tree_n,
  #                          trControl = train.control, ...)
  node.mod <- caret::train(x = node.Data,
                           y = node.ClassLabels,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = TRUE,
                           proximity = FALSE,
                           outscale = FALSE,
                           preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)


  return(node.mod)
}


NodeTrainer4 <- function(Rdata, tree, node, f_n=200, tree_n=500, PC_n=40, ...){

  node.Data <- SubsetTData(Tdata = Rdata, tree = tree, node = node)
  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  #node.Data <- node.Data[, apply(node.Data, 2, var) != 0]

  node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)
  if(labs_l[node] != "TaxaRoot"){
    childNodes <- GetChildNodeLeafs(tree = tree, node = node)
    childLeafs <- NULL
    for (i in which(lapply(childNodes, length) > 0)){
      childLeafs <- c(childLeafs, childNodes[i][[1]])
    }
    outGroupLeafs <- tree$tip.label[!tree$tip.label %in% childLeafs]
    SampSize.min <- min(table(Rdata$ClassLabels)[outGroupLeafs])
    for(x in outGroupLeafs){
      node.outData <- droplevels(Rdata[which(Rdata$ClassLabels == x), ])
      node.outData <- node.outData[sample(rownames(node.outData), size = SampSize.min), ]
      node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
      node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
      node.Data <- rbind(node.Data, node.outData)
    }
  }else{
    RandN <- median(table(node.ClassLabels))
    node.outData <- droplevels(Rdata[sample(rownames(Rdata), size=RandN), ])
    node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
    node.ClassLabels <- node.Data$ClassLabels
    node.outData <- droplevels(subset(node.outData, select=-c(ClassLabels)))
    node.outData <- Shuffler(df = node.outData)
    node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
    node.Data <- rbind(node.Data, node.outData)
  }

  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))

  #Here transform the node.Data into PC space:
  #node.Data <- node.Data[, apply(node.Data, 2, var) != 0]

  #Select top highly varible genes:
  vmr <- Seurat::FindVariableFeatures(t(expm1(node.Data)),verbose=F)
  vmr <- vmr[order(vmr$vst.variance.standardized, decreasing = T),]
  node.Data <- node.Data[, head(rownames(vmr), 2000)]

  pcatrain <- prcomp(node.Data, center = TRUE, scale=TRUE, rank. = PC_n)
  PCs.sig <- selectSigPCs(pcatrain, node.ClassLabels)

  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)

  node.mod <- caret::train(x = pcatrain$x[,PCs.sig],
                           y = node.ClassLabels,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = TRUE,
                           proximity = FALSE,
                           outscale = FALSE,
                           #preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)

  node.mod[["trainPCA"]] <- pcatrain
  return(node.mod)
}

NodeTrainer5 <- function(Rdata, tree, node, f_n=200, tree_n=500, ...){

  node.Data <- SubsetTData(Tdata = Rdata, tree = tree, node = node)
  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  #node.Data <- node.Data[, apply(node.Data, 2, var) != 0]

  node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)
  if(labs_l[node] != "TaxaRoot"){
    childNodes <- GetChildNodeLeafs(tree = tree, node = node)
    childLeafs <- NULL
    for (i in which(lapply(childNodes, length) > 0)){
      childLeafs <- c(childLeafs, childNodes[i][[1]])
    }
    outGroupLeafs <- tree$tip.label[!tree$tip.label %in% childLeafs]
    SampSize.min <- min(table(Rdata$ClassLabels)[outGroupLeafs])
    for(x in outGroupLeafs){
      print(x)
      node.outData <- droplevels(Rdata[which(Rdata$ClassLabels == x), ])
      node.outData <- node.outData[sample(rownames(node.outData), size = SampSize.min), ]
      node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
      node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
      node.Data <- rbind(node.Data, node.outData)
    }
  }


  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))

  # Pool <- NULL
  # for(cc in names(tabl)){
  #   #SampSize.min instead of tabl[cc]
  #   nonzeroFeatures <- names(node.Data)[colSums(node.Data[node.ClassLabels == cc, ] != 0)/SampSize.min > .5]
  #   Pool <- append(Pool, nonzeroFeatures)
  # }
  # Pool <- unique(sort(Pool))
  #
  # node.Data <- droplevels(subset(node.Data, select=c(Pool)))


  #Select top highly varible genes:
  vmr <- Seurat::FindVariableFeatures(t(expm1(node.Data)),verbose=F)
  vmr <- vmr[order(vmr$vst.variance.standardized, decreasing = T),]
  P_dict <- head(rownames(vmr), f_n)

  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))


  #node.Data$ClassLabels <- node.ClassLabels
  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE,
                                       sampling = "down",
                                       preProcOptions = list(thresh = 0.95,
                                                             ICAcomp = 3,
                                                             k = 5,
                                                             freqCut = 95/5,
                                                             uniqueCut = 10,
                                                             cutoff = 0.9))

  node.mod <- caret::train(x = node.Data,
                           y = node.ClassLabels,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = TRUE,
                           proximity = FALSE,
                           outscale = FALSE,
                           preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)


  return(node.mod)
}


#' @param
#' @usage  genePairs <- FindPairs(P_dict, node.Data, node.ClassLabels)
PairsPool <- function(P_dict, node.Data, node.ClassLabels, score_threshold=0.5){
  library(tspair)
  genes <- P_dict
  df <- t(node.Data)

  genePairs <- c()
  pairs.pool <- c()

  for(c in unique(as.character(node.ClassLabels))){
    genes <- P_dict
    df <- t(node.Data)
    #create the gene pairs:
    data <- df[P_dict, ]

    node.labs <- as.character(node.ClassLabels)
    node.labs[which(node.labs != c)] <- "A"
    node.labs[which(node.labs == c)] <- "B"

    while(length(genes) >= 2){
      tsp1 <- tspcalc(as.matrix(data), node.labs)
      pairs <- rownames(tsp1$tspdat)
      if((!pairs[1] == pairs[2]) && (tsp1$tspscore > score_threshold)){
        pairs.pool <- append(pairs.pool,pairs[1])
        pairs.pool <- append(pairs.pool,pairs[2])
      }
      genes <- genes[!genes %in% pairs]
      data <- data[genes,]
    }
  }
  return(unique(pairs.pool))
}



#' A function to tally votes from the model locally trained.
#' @param model
#' @param QueData is the prepared data matrix ready to be used in predict
#' @param format type of prediction output, "prob" or "resp".
CalibratedPvoteR_caret <- function(model, mlr, QueData, format="prob", node=NULL){
  if(is.null(node)){
    QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
  } else{
    if(format == "prob"){
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
      colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
    } else {
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "raw", scale=T, center=T))
      colnames(QuePvotes) <- as.character(node)
    }
  }
  print(head(QuePvotes))

  #Calibrate the probabilities:
  QuePvotes.cal <- predict(mlr, newdata=QuePvotes, type="prob")
  print(head(QuePvotes.cal))

  if(dim(QuePvotes)[2] < 3){
    print("dim(QuePvotes)[2] < 3")
    QuePvotes.cal <- data.frame(QuePvotes.cal, 1-QuePvotes.cal)
    colnames(QuePvotes.cal) <- colnames(QuePvotes)#Double-check the colnames.
  }
  Calibrated_and_Scaled <- data.frame(row.names = rownames(QuePvotes.cal))
  for(class in colnames(QuePvotes.cal)){
    vec <- QuePvotes.cal[, class]
    vec <- (vec - min(vec))/(max(vec) - min(vec))
    Calibrated_and_Scaled <- cbind(Calibrated_and_Scaled, data.frame(class=vec))
  }
  colnames(Calibrated_and_Scaled) <- colnames(QuePvotes.cal)
  print(head(Calibrated_and_Scaled))
  return(Calibrated_and_Scaled)
}


GetSigMods <- function(RefData, ClassLabels, refmod, ...){

  Cal.data <- NoiseInject(RefData = RefData, ClassLabels = ClassLabels, refmod = refmod)
  Pvotes <- ProbTraverser(Query = Cal.data[["data"]], refmod = refmod)

  votes <- data.frame(Prior=Cal.data[["Y"]], Pvotes)
  classes <- unique(votes$Prior)
  labs_l <- c(refmod@tree[[1]]$tip.label, refmod@tree[[1]]$node.label)
  sigmoidMods <- list()

  for(n in 1:length(labs_l)){
    if(labs_l[n] == "TaxaRoot"){next}
    Leafs <- NULL
    if(labs_l[n] %in% refmod@tree[[1]]$tip.label){
      Leafs <- labs_l[n]
    }else{
      leaf <- GetChildNodeLeafs(tree = refmod@tree[[1]], node = n)
      leaf <- leaf[which(lapply(leaf, length)>0)]
      for(j in 1:length(leaf)){
        Leafs <- append(Leafs, leaf[[j]])
      }
    }
    votes <- votes[order(votes[[labs_l[n]]]),]
    Labs <- ifelse(votes$Prior %in% Leafs, 1, 0)
    dfr <- data.frame(votes[[labs_l[n]]], Labs)
    colnames(dfr) <- c("x","y")
    # training a logistic regression model on the cross validation dataset
    glmmod <- glm(y~x, data = dfr, family = binomial, maxit = 100, y=FALSE, model=FALSE)
    sigmoidMods[[labs_l[n]]] <- stripGlmModel(glmmod)
  }
  outs <- grep("OutGroup", names(votes), value = T)
  for(x in 1:length(outs)){
    int.node <- gsub("_OutGroup", "", outs[x])
    n <- match(int.node, labs_l)
    Leafs <- NULL
    leaf <- GetChildNodeLeafs(tree = refmod@tree[[1]], node = n)
    leaf <- leaf[which(lapply(leaf, length)>0)]
    for(j in 1:length(leaf)){
      Leafs <- append(Leafs, leaf[[j]])
    }
    votes <- votes[order(votes[[outs[x]]]),]
    Labs <- ifelse(votes$Prior %in% Leafs, 0, 1)#Labels (1,0) reversed since outgroup
    dfr <- data.frame(votes[[outs[x]]], Labs)
    colnames(dfr) <- c("x","y")
    # training a logistic regression model on the cross validation dataset
    glmmod <- glm(y~x, data = dfr, family = binomial, maxit = 100, y=FALSE, model=FALSE)
    sigmoidMods[[outs[x]]] <- stripGlmModel(glmmod)
  }

  return(sigmoidMods)
}



DetermineAlpha <- function(refMod, RefData, Prior){
  refMod@alpha[[1]] <- 0
  Hierobj <- HieRFIT(Query = RefData, refMod = refMod, Prior = Prior)
  fails <- Hierobj@Evaluation$Projection != Hierobj@Prior
  alpha <- mean(Hierobj@Evaluation[fails,]$Certainty)
  if(is.na(alpha)){alpha<-0}
  return(alpha)
}

DetermineAlpha3 <- function(refMod, RefData, Prior){
  tree <- refMod@tree[[1]]
  labs_l <- c(tree$tip.label, tree$node.label)
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"]

  alpha.list.def <- list()
  for(node.lab in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    alp <- rep(0, length(AncPath))
    names(alp) <- AncPath
    alpha.list.def[[node.lab]] <- alp
  }

  refMod@alpha[[1]] <- alpha.list.def

  Hierobj <- HieRFIT(Query = RefData, refMod = refMod, Prior = Prior)

  succs <- Hierobj@Evaluation$Projection == Hierobj@Prior
  succsLabs <- as.character(Hierobj@Evaluation[succs,]$Projection)
  succsCerst <- Hierobj@CertaintyValues[succs,]

  fails <- Hierobj@Evaluation$Projection != Hierobj@Prior
  failsLabs <- as.character(Hierobj@Evaluation[fails,]$Projection)
  failsCerst <- Hierobj@CertaintyValues[fails,]

  alpha.list <- list()
  for(node.lab in labs_l[labs_l %in% failsLabs]){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    al.list <- NULL
    for(i in 1:length(AncPath)){
      if(any(failsLabs == node.lab)){
        maxcer <- mean(failsCerst[failsLabs == node.lab, AncPath[i]])
      }else{
        maxcer <- 0
      }
      al.list <- append(al.list, maxcer)
    }
    names(al.list) <- AncPath
    alpha.list[[node.lab]] <- al.list
  }

  for(node.lab in labs_l[!labs_l %in% failsLabs]){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    al.list <- NULL
    for(n in 1:length(AncPath)){
      list.n <- NULL
      for(i in 1:length(alpha.list)){
        if(AncPath[n] %in% names(alpha.list[[i]])){
          n.cer <- alpha.list[[i]][AncPath[n]]
        }else{
          n.cer <- 0
          names(n.cer) <- AncPath[n]
        }
        list.n <- append(list.n, n.cer)
      }
      list.n <- mean(list.n, na.rm = T)
      al.list <- append(al.list, list.n)
    }
    names(al.list) <- AncPath
    alpha.list[[node.lab]] <- al.list
  }

  return(alpha.list)
}



#' A function to evaluate performance of HieRFIT with various Certainty thresholds.
#' @param RefSeuObj reference seurat object.
#' @param IdentityCol the column 'number' in the metadata slot showing the cell type identities.
#' @param refMod
#' @param Uinter number of intervals to assess certainty threshold.
#' @param perm_n # of permutations.
#' @param samp_n # samples to draw at each permutation.
EvaluateCertainty <- function(RefSeuObj, IdentityCol, refMod, Uinter=20, perm_n=10, samp_n=100){
  hPRF_table <- numeric()
  for(p in 1:perm_n){
    testcells <- sample(rownames(RefSeuObj@meta.data), samp_n)
    for( ai in seq(1:Uinter)/Uinter){
      testObj <- HieRFIT(Query = as.matrix(RefSeuObj@data)[, testcells],
                         refMod = refMod,
                         Prior = RefSeuObj@meta.data[testcells, IdentityCol],
                         alpha = ai)

      PriorPostTable <- data.frame(Prior = testObj@Prior,
                                   Projection = testObj@Evaluation$Projection)
      hPRF.out <- hPRF(tpT = PriorPostTable, tree = refMod@tree[[1]])
      hPRF.out <- c(Perm=p,
                    U_threshold=ai,
                    hPRF.out,
                    Undetermined=nrow(PriorPostTable[PriorPostTable$Projection == "Undetermined",])/length(PriorPostTable[,1]))
      print(hPRF.out)
      hPRF_table <- rbind(hPRF_table, hPRF.out)
    }
  }
  return(hPRF_table)
}



if(package == 'caret'){
  train.control <- caret::trainControl(method="oob",
                                         returnData = FALSE,
                                         savePredictions = "none",
                                         returnResamp = "none",
                                         allowParallel = TRUE,
                                         classProbs =  TRUE,
                                         trim = TRUE)
    print("training with caret...")
    node.mod <- caret::train(x = node.Data,
                             y = node.ClassLabels,
                             method = "rf",
                             norm.votes = TRUE,
                             importance = TRUE,
                             proximity = FALSE,
                             outscale = FALSE,
                             preProcess = c("center", "scale"),
                             ntree = tree_n,
                             trControl = train.control, ...)
  }else{
    print("training with scikit-learn...")
    f_sub <- reduceFeatures(x = node.Data,
                            y = node.ClassLabels,
                            ntree = tree_n)
    node.Data <- droplevels(subset(node.Data, select=c(f_sub)))

    node.mod <- ski_train(x=node.Data,
                          y=node.ClassLabels,
                          ntree = tree_n)
  }



reduceFeatures <- function(x, y, ntree=500, cutoff=0.8, plot=FALSE){

  node.mod <- ski_train(x=x,
                        y=y,
                        ntree = ntree)

  feats <- node.mod$feature_importances_
  names(feats) <- colnames(x)
  feats <- sort(feats, decreasing = T)
  if(plot == TRUE){
    plt <- plot(cumsum(feats), pch=20, xlab = "features")
    abline(h=0.8, col="red", lwd=3, lty=2)
  }
  feats <- names(feats[cumsum(feats) < 0.8])

  return(feats)
}

ski_train <- function(x, y, ntree){
  #reticulate::use_python(python = "/Users/yasinkaymaz/miniconda3/bin/python")
  mod <- skiRFC_trainer(x=x,
                        y=y,
                        tree_n=as.integer(ntree) )
  f_set <- mod$feature_importances_
  names(f_set) <- colnames(x)
  f_set <- sort(f_set, decreasing = T)
  #names(head(feat1, 1000))
  mod$finalModel <- NULL
  mod$finalModel$xNames <- colnames(x)
  mod$finalModel$classes <- mod$classes_
  return(mod)
}



NoiseInject <- function(RefData, ClassLabels, refMod=NULL){
  NoisedRef <- list()
  Caldata <- RefData
  rownames(Caldata) <- FixLab(rownames(Caldata))#Not necessary
  features <- NULL
  if(is.null(refMod)){
    features <- rownames(Caldata)
  }else{
    for( i in names(refMod@model)){ features <- append(features, refMod@model[[i]]$finalModel$xNames)}
  }
  features <- unique(features)
  Caldata <- Caldata[features, ]
  CalY <- FixLab(ClassLabels)
  Caldata.noised <- Caldata

  # max(table(CalY))  #For balancing the class data
  # Caldata[, CalY == unique(CalY)[1]]

  rands <- seq(1:5)
  CaldataFull <- Caldata
  for(r in rands){
    for(x in 1:length(Caldata[1,])){
      nz_idx <- which(!Caldata[,x] %in% c(0))
      if(r > length(nz_idx)){
        r = length(nz_idx)-1
      }
      idx.to.set.zero <- sample(nz_idx, r)
      Caldata.noised[idx.to.set.zero, x] <- 0
    }
    CaldataFull <- cbind(CaldataFull, Caldata.noised)
  }

  TaxaOutdata <- Caldata[, sample(colnames(Caldata), size = 500)]#Replace 500 with an adaptive.
  TaxaOutdata <- Shuffler(df = TaxaOutdata)
  CaldataFull <- cbind(CaldataFull, TaxaOutdata)

  colnames(CaldataFull) <- make.names(colnames(CaldataFull), unique = TRUE)
  CalYFull <- c(rep(CalY, length(rands)+1), rep("TaxaOut", 500))

  NoisedRef[["data"]] <- CaldataFull
  NoisedRef[["Y"]] <- CalYFull

  return(NoisedRef)
}


selectSigPCs <- function(pcatrain, ClassLabels){
  pcadata <- data.frame(pcatrain$x, ClassLabels = ClassLabels)
  cls <- levels(ClassLabels)
  PC_n <- length(pcatrain$x[1,])
  ptab <- NULL
  for(i in 1:PC_n){
    PC.stats <- NULL
    for(c in cls){
      p.cl <- t.test(pcadata[pcadata$ClassLabels == c, i],
                     pcadata[pcadata$ClassLabels != c, i])$p.value
      PC.stats <- c(PC.stats, p.cl)
    }
    names(PC.stats) <- cls
    pc.col <- paste("PC", i, sep = "")
    ptab <- cbind(ptab, pc.col=PC.stats)
  }

  ptab <- as.data.frame(ptab)
  colnames(ptab) <- colnames(pcadata)[-length(pcadata[1,])]
  ptab <- ptab*length(cls)*PC_n#Correct for the multiple test. Bonferroni.
  #Select only PCs which are significanly separating at least one of the class.
  PCs.sig <- colnames(ptab[, apply(ptab < 0.05, 2 ,any)])

  if(length(PCs.sig) == 0){PCs.sig <- paste(rep("PC", 3), 1:3, sep="")}
  return(PCs.sig)
}



#' A function to calculate class scores for all internal and tip node classes.
#' @param tree
#' @param nodes_P_all a table output from CTTraverser() which contains all class probabilities from every node.
#' @return P_path_prod a table for products of ancestor node scores.
ClassProbCalculator <- function(tree, nodes_P_all){
  P_path_prod <- data.frame(row.names = rownames(nodes_P_all))
  clabs <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  for(cl in clabs){
    map <- GetAncestPath(tree = tree, class = cl)
    if(length(map) > 1){
      nt_prob <- data.frame(apply(nodes_P_all[, map], 1, prod))
    }else{
      nt_prob <- data.frame(nodes_P_all[,map])
    }
    colnames(nt_prob) <- cl
    P_path_prod <- cbind(P_path_prod, nt_prob)
  }
  return(P_path_prod)
}


ski_predict <- function(fit, newdata, type='prob'){
  #reticulate::use_python(python = "/Users/yasinkaymaz/miniconda3/bin/python")
  #Sys.which("python")
  probs <- skiRFC_predict(fit=fit, x=newdata, outcome=type)
  if(type == 'prob'){
    colnames(probs) <- fit$classes_
    rownames(probs) <- rownames(newdata)
  }
  return(probs)
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  #This function has been copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


tree2table <- function(refMod){
  allNodeNames <- c(refMod@tree[[1]]$tip.label,refMod@tree[[1]]$node.label)
  treeTable<-cbind(allNodeNames[refMod@tree[[1]]$edge[,1]],allNodeNames[refMod@tree[[1]]$edge[,2]])
  rootNode<-treeTable[1,1]
  treeTable<-t(apply(treeTable,1,function(row){
    newRow<-setdiff(row,rootNode)
    newRow<-gsub("_"," ",newRow)
    c(newRow,rep("",length(row)-length(newRow)))
  }))
  treeTable<-as.data.frame(treeTable)
  return(treeTable)
}



PlotPredictions <- function(SeuratObject, model, save.pdf=T, outputFilename="plotpredictions") {
  #Evaluate model prediction accuracy:
  conf.mat <- model$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>%
    mutate(freq = 100*value/sum(value))

  class.n = length(model$classes)

  pdf(paste(outputFilename,".pdf",sep=""),width= 1.5*class.n, height = 1.5*class.n)


  require(gridExtra)

  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +
    geom_tile(color = "white")+
    coord_equal()+
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

PlotPredictions2 <- function(SeuratObject, model, priorLabels, outputFilename="plotpredictions") {
  #Evaluate model prediction accuracy:
  conf.mat <- model$finalModel$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>% mutate(freq = 100*value/sum(value))
  class.n = length(model$finalModel$classes)
  errorSize <- as.data.frame(cbind(model$finalModel$confusion[,"class.error"], head(colSums(model$finalModel$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- getTrainPerf(model)["TrainAccuracy"]*100
  di <- round(sqrt(class.n),digits = 0)+1

  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +
    geom_tile(color = "white")+
    coord_equal()+
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


PlotPredictions22 <- function(SeuratObject, outputFilename="plotpredictions") {

  #Prediction Performance
  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=FinalPrediction,fill=FinalPrediction),stat = "count")+
    geom_violin(aes(x=FinalPrediction,y=FinalBestProb*max(table(SeuratObject@meta.data$FinalPrediction)),fill=FinalPrediction))+
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$FinalPrediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))

  p4 <- TSNEPlot(SeuratObject, group.by="FinalPrediction",do.label=T, do.return = T)

  p3p4 <- cowplot::plot_grid(p3,p4, ncol = 1, nrow=2)
  save_plot(filename = paste(outputFilename,".prediction-stats.pdf",sep=""),plot = p3p4,base_height = 20, base_width = 24)

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




PlotProjectionStats <- function(SeuratObject, outputFilename="plotpredictions") {
  pdf(paste(outputFilename,".pdf",sep=""), width= 12, height = 8)

  TSNEPlot(SeuratObject)

  TSNEPlot(SeuratObject, group="Projection", do.label = T)

  p <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=Projection,fill=Projection),stat = "count")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Projection outcome", sep=""))
  print(p)

  p <- ggplot(data=SeuratObject@meta.data)+
    geom_boxplot(aes(x=Projection,y=Score,fill=Projection))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Classification Score",title=paste("Prediction outcome", sep=""))
  print(p)

  p <- ggplot(data=SeuratObject@meta.data)+
    geom_boxplot(aes(x=Projection,y=Certainty,fill=Projection))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Classification Certainty",title=paste("Prediction outcome", sep=""))
  print(p)

  dev.off()
}



#' This function calculates a scaled Kullback-Leibler divergence
#' @param probs a list of observed probability scores.
#' @return KLscaled = KLe/KLmax, where KLe is the empirical divergence given
#' the distributions while KLmax is the maximum KLe value that can be achieved
#' given the number of items.
KLeCalc <- function(probs){
  class_n <- length(probs)
  null <- c(1, rep(0, class_n-1))
  KLe <- entropy::KL.empirical(y1 = as.numeric(probs), y2 = rep(1/class_n, class_n))
  KLmax <- entropy::KL.empirical(y1 = as.numeric(null), y2 = rep(1/class_n, class_n))
  KLscaled <- KLe/KLmax
  return(KLscaled)
}


#' A function to calculate Asymmetric Entropy
#' @param p measured probability array
#' @param w empirically calculated W set
#' @return U set of certainty values for each probility outcome in p given w.
GetCertaintyArray <- function(p, w){
  U <- numeric()
  for(i in 1:length(p)){
    w[i] <- w[i]+1e-10#to prevent math err.
    if(p[i] > w[i]){
      lambda=1
    }else{
      lambda=-1
    }
    U <- c(U, (lambda*(p[i]-w[i])^2)/( ((1-2*w[i])*p[i])+w[i]^2 ) )
  }
  U <- as.numeric(U)
  return(U)
}


TipBias <- function(tree, confmat){
  #confmat is table() of Prior/Prediction comparison
  err.rate <- NULL
  for(i in 1:length(tree$tip.label)){
    nn <- length(GetAncestorsPath(tree=tree, i)[[2]])
    leafErr <- 1-confmat[tree$tip.label[i],tree$tip.label[i]]/sum(confmat[tree$tip.label[i],])
    print(paste(i,nn,leafErr,sep = "    "))
    err.rate <- c(err.rate, leafErr)
  }
  err.rate <- data.frame(err.rate)
  rownames(err.rate) <- tree$tip.label
  return(err.rate)
}




#' An internal function to prepare a training/test dataset for model generation.
#' @param Data a Normalized expression data matrix, genes in rows and samples in columns.
#' @param Predictors the predictor feature list selected by FeatureSelector.
#' @param ClassLabels [optional] A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param alpa variation cutoff for filtering data
#' @keywords data preparation
#' @export
#' @usage trainingData <- DataReshaper(Data = as.matrix(SeuratObject@data), Predictors = genes, ClassLabels = SeuratObject@meta.data$CellTypes)
DataReshaper <- function(Data, Predictors, ClassLabels, alpa=0.1, ...) {
  if(missing(Predictors)){#prepare the data for PCA
    TData <- as.data.frame(t(Data))
    indx <- sapply(TData, is.factor)
    TData[indx] <- lapply(TData[indx], function(x) as.numeric(as.character(x)))
    #Filter candidate predictors before PCA
    TData <- TData[, apply(TData, 2, var) != 0]
    TData <- droplevels(TData[, which(matrixStats::colVars(as.matrix(TData)) > alpa)])
    return(TData)

  } else {

    if (missing(ClassLabels)) {
      #Then the output is for Query
      QueData <- as.data.frame(t(Data), col.names=rownames(Data))
      colnames(QueData) <- make.names(colnames(QueData))
      QueData_sub <- droplevels(QueData[, which(colnames(QueData) %in% Predictors)])
      #missing Predictors
      mp <- Predictors[which(!Predictors %in% colnames(QueData))]
      #Add missing predictors into QueData by setting to 0.
      mp_df <- data.frame(matrix(0,
                                 ncol = length(mp),
                                 nrow = length(colnames(Data))))
      colnames(mp_df) <- mp
      QueData <- cbind(QueData_sub, mp_df)
      QueData <- QueData[, Predictors]
      return(QueData)

    } else {

      RefData <- as.data.frame(t(Data), col.names=rownames(Data))
      colnames(RefData) <- make.names(colnames(RefData))
      #convert factors to numeric
      indx <- sapply(RefData, is.factor)
      RefData[indx] <- lapply(RefData[indx], function(x) as.numeric(as.character(x)))
      RefData$ClassLabels <- factor(make.names(ClassLabels))
      RefData <- droplevels(RefData[, c(Predictors, "ClassLabels")])

      return(RefData)
    }
  }
}

#' A wrapper function for random forest.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(RefData).
#' @param mod.meth The model training method, "rf" for random forest.
RandForestWrap <- function(RefData, ClassLabels, prefix, mod.meth, train.control, thread=5, ...){
  #library(doParallel)
  library(caret)
  library(doMC)
  registerDoMC(cores = thread)
  #1. Select the predictors.
  P_dicts <- FeatureSelector(Data = RefData,
                             ClassLabels = ClassLabels,
                             PCs = 10,
                             num = 200,
                             prefix = prefix,
                             doPlots = F, ...)
  #2. Prepare the reference data.
  TrainData <- DataReshaper(Data = RefData,
                            Predictors = P_dicts,
                            ClassLabels = ClassLabels, ...)
  print(RefData[1:5, 1:5])
  #3. Train the model.
  #cl <- makePSOCKcluster(ncor)
  #registerDoParallel(cl)
  model <- caret::train(ClassLabels~.,
                        data = TrainData,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = FALSE,
                        proximity = FALSE,
                        outscale = FALSE,
                        preProcess = c("center", "scale"),
                        ntree=50,
                        trControl = train.control)
  #stopCluster(cl)

  return(model)
}

#' A wrapper function for random forest.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param mod.meth The model training method, "svmLinear" for support vector machine.
SvmWrap <- function(RefData, ClassLabels, prefix, mod.meth, train.control, ncor=5, ...){
  library(doParallel)
  #1. Select the predictors.
  P_dicts <- FeatureSelector(Data = RefData,
                             ClassLabels = ClassLabels,
                             PCs = 10,
                             num = 200,
                             prefix = prefix,
                             doPlots = F, ...)
  #2. Prepare the reference data.
  TrainData <- DataReshaper(Data = RefData,
                            Predictors = P_dicts,
                            ClassLabels = ClassLabels, ...)

  #3. Train the model.
  cl <- makePSOCKcluster(ncor)
  registerDoParallel(cl)
  model <- caret::train(ClassLabels~.,
                        data = TrainData,
                        trControl = train.control,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = TRUE,
                        proximity = TRUE,
                        preProcess = c("center", "scale"),
                        tuneLength = 10, ...)
  stopCluster(cl)

  return(model)
}

#' A wrapper function for quick data processing with Seurat functions
#' This function allows to determine highly variable genes and scale their expression,
#' run PCA, tSNE, and cluster detection.
#' @param SeuratObj a Seurat S4 object.
#' @param scale.only.var a boolean variable to determine whether to scale entire data or only highly variable genes. Default is True.
#' @param PCs number of PCs to be included in the downstream analysis
#' @param vars2reg variables to be regressed out when scaling the expression data.
#' @param perp perplexity parameter to be passed to RunTSNE function from Seurat.
#' @keywords seurat quick
#' @export
#' @usage pbmc <- QuickSeurat(pbmc, scale.only.var=F, PCs=5, perp=20)
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

    #Save seurat object:
    save(SeuratObj, file=paste(ProjectLabel, ".seurat.Robj", sep=""))
  }

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

FlatRF <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions") {

  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
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

  if (!is.nan(mmDGm)) {
    if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
      warning("A significant portion of features are missing...")
    }
  }


  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class.n <- length(model$finalModel$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="raw")
  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

  if (missing(priorLabels)) {
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
      scale_x_discrete(limits = c("Prior", "Clusters", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
  }

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function


#' An internal function to collect model training parameters and direct them to model creation.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param mod.meth The model training for hierarchical random forest. Default is "hrf"
#' @param cv.k Fold cross validation. Default is 5.
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param save.int.f Boolean to save model. Default is False.
#' @keywords
#' @export
#' @usage rf.model <- Modeller(RefData = as.matrix(pbmc1@data), ClassLabels = pbmc1@meta.data$ClusterNames_0.6)
Modeller <- function(RefData, ClassLabels=NULL, mod.meth="rf", thread=NULL, tree=NULL, save.int.f=FALSE, ...){
  if(mod.meth == "hrf"){
    try(if(missing(tree)|| missing(ClassLabels) || missing(RefData))
      stop("Please, provide the required inputs!"))
    model <- HieRandForest(RefData = RefData,
                           ClassLabels = ClassLabels,
                           tree, thread = thread)
  }

  return(model)
}


#' Homology mapping via orthologous genes between mouse and rat.DEPRICATED
Gmor <- function(RatGenes){
  # This function retrieves mouse homolog associated gene names of Rat genes.
  #library(biomaRt)
  ensembl.rat <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl",
                                  host = "www.ensembl.org",
                                  ensemblRedirect = FALSE)
  R2M.ort <- biomaRt::getBM(attributes = c("external_gene_name",
                                  "mmusculus_homolog_associated_gene_name",
                                  "mmusculus_homolog_orthology_confidence",
                                  "mmusculus_homolog_orthology_type",
                                  "ensembl_gene_id",
                                  "mmusculus_homolog_ensembl_gene"),
                   filters = 'external_gene_name',
                   values = RatGenes,
                   uniqueRows = T,
                   mart = ensembl.rat)
  R2M.ort <- R2M.ort[which(R2M.ort$mmusculus_homolog_orthology_confidence == "1"), ]
  R2M.ort <- R2M.ort[which(R2M.ort$mmusculus_homolog_orthology_type == "ortholog_one2one"), ]
  return(R2M.ort)
}

RAMmerger <- function(RatObj, MouseObj){
  ort <- Gmor(RatGenes = rownames(RatObj@data))
  mm.data <- as.matrix(MouseObj@raw.data[which(rownames(MouseObj@raw.data) %in% ort$mmusculus_homolog_associated_gene_name),])
  rownames(mm.data) <- ort[match(rownames(mm.data), ort$mmusculus_homolog_associated_gene_name),]$external_gene_name
  rownames(mm.data) <- make.names(rownames(mm.data), unique = T)
  Mouse <- SeuratWrapper(ExpData = mm.data, ProjectLabel = "Mouse_data", NewMeta = MouseObj@meta.data, Normalize = T, dump.files = F)

  Rat.data <- as.matrix(RatObj@raw.data[ort[which(ort$mmusculus_homolog_associated_gene_name %in% rownames(zeisub.data)),]$external_gene_name,])
  rownames(Rat.data) <- make.names(rownames(Rat.data), unique = T)
  Rat <- SeuratWrapper(ExpData = Rat.data, ProjectLabel = "Rat_data",  NewMeta = RatObj@meta.data, Normalize = T, dump.files = F)

  ccaMergedMouseRat <- SeuratCCAmerger(listofObjects = c(Mouse, Rat))

  return(ccaMergedMouseRat)
}



stripGlmModel = function(cm) {
  cm$y = c()
  cm$model = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  return(cm)
}


CVRunner <- function(Ref, ClassLabels, TreeTable=NULL, cv_k=5, method="hrf"){

  # Create K-fold data split:
  CMtab <- vector("list", length = cv_k)
  flds <- caret::createFolds(ClassLabels, k = cv_k, list = TRUE, returnTrain = FALSE)
  hPRF_cv <- do.call(rbind.data.frame,
                     lapply(1:cv_k,
                            function(i){
                              #Create hieR
                              trainClassLabels <- ClassLabels[-flds[[i]]]
                              trainRef <- Ref[, -flds[[i]]]

                              refmod <- CreateHieR(RefData = trainRef,
                                                   ClassLabels = trainClassLabels,
                                                   TreeTable = TreeTable)
                              #Hierfit
                              testClassLabels <- ClassLabels[flds[[i]]]
                              testRef <- Ref[, flds[[i]]]


                              testObj <- HieRFIT(Query = testRef, refMod = refmod, Prior = testClassLabels)
                              PriorPostTable <- data.frame(Prior=testObj@Prior, Projection = testObj@Evaluation$Projection)
                              hPRF.out <- hPRF(tpT = PriorPostTable, tree = refmod@tree[[1]])
                              hPRF.out <- c(Tool="HieRFIT",
                                            CV_k=i,
                                            hPRF.out)

                              print(hPRF.out)
                              cc <- t(hPRF.out)
                              return(cc)
                            }
                     )
  )
  return(hPRF_cv)
}

Rander <- function(df){
  if(dim(df)[1] > 1000){Sn <- 1000}else{Sn <- dim(df)[1]}
  dfr <- as.matrix(df[sample(row.names(df), size = Sn),])
  for(i in 1:length(dfr[1,])){#for each col:
    for(j in 1:length(dfr[,1])){#for each row:
      dfr[,i] <- sample(dfr[,i])
      dfr[j,] <- sample(dfr[j,])
    }
  }
  return(as.data.frame(dfr))
}

RandomizeR <- function(df, n=10){
  #set.seed(192939)
  dfRand <- NULL
  for (i in 1:n){
    dfR <- df[sample(nrow(df)), sample(ncol(df))]
    rownames(dfR) <- rownames(df)
    colnames(dfR) <- colnames(df)
    dfRand <- rbind(dfRand, dfR)
  }
  dfRand <- dfRand[sample(rownames(dfRand),size = nrow(df)),]
  return(dfRand)
}

FakeRandomizeR <- function(df, seed.num=192939){
  set.seed(seed.num)
  dfRand <- df
  rownames(dfRand) <- sample(rownames(df))
  colnames(dfRand) <- sample(colnames(df))
  return(dfRand)
}



evaluate <- function(TrueLabels, PredLabels, Indices = NULL, HierModPath=NULL){
  "
  This function was taken from https://github.com/tabdelaal/scRNAseq_Benchmark
  "
  #true_lab <- unlist(read.csv(TrueLabelsPath))
  #pred_lab <- unlist(read.csv(PredLabelsPath))
  true_lab <- TrueLabels
  pred_lab <- PredLabels

  true_lab <- FixLab(xstring = true_lab)
  pred_lab <- FixLab(xstring = pred_lab)

  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }

  if(!is.null(HierModPath)){

    if(class(HierModPath)[1] == "RefMod"){
      refmod <- HierModPath
    }else{
      suppressPackageStartupMessages(library(HieRFIT))
      refmod <- LoadHieRMod(fileName=HierModPath)
    }

    hPRFtab <- hPRF(tpT = as.data.frame(cbind(true_lab, pred_lab)), tree = refmod@tree[[1]])

    }else{hPRFtab <- NULL}

  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))

  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)

  pred_lab = gsub('Node..','Node',pred_lab)

  conf_F1 <- table(true_lab,pred_lab,exclude = c('Undetermined',
                                                 'unassigned',
                                                 'Unassigned',
                                                 'Unknown',
                                                 'rand',
                                                 'Node',
                                                 'Int.Node',
                                                 'ambiguous',
                                                 'unknown'))

  F1 <- vector()
  sum_acc <- 0

  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
  }

  pop_size <- pop_size[pop_size > 0]

  names(F1) <- names(pop_size)

  med_F1 <- median(F1)
  mean_F1 <- mean(F1)

  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == 'Undetermined') +
    sum(pred_lab == 'unassigned') +
    sum(pred_lab == 'Unassigned') +
    sum(pred_lab == 'rand') +
    sum(pred_lab == 'Unknown') +
    sum(pred_lab == 'unknown') +
    sum(pred_lab == 'Node') +
    sum(pred_lab == 'Int.Node') +
    sum(pred_lab == 'ambiguous')
  per_unlab <- num_unlab / total

  num_Interlab <- sum(pred_lab == 'Node') +
    sum(pred_lab == 'Int.Node')
  per_Interlab <- num_Interlab / total

  acc <- sum_acc/sum(conf_F1)

  result <- list(Conf = conf, MeanF1=mean_F1, MedF1 = med_F1, F1 = F1, Acc = acc, PercInter= per_Interlab, PercUnl = per_unlab, PopSize = pop_size, hPRF=hPRFtab)

  return(result)
}


NodeTrainer2 <- function(Rdata, tree, node, f_n=200, tree_n=500, switchBox='off', ...){

  node.Data <- SubsetTData(Tdata = Rdata, tree = tree, node = node)
  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  node.Data <- node.Data[, apply(node.Data, 2, var) != 0]

  if(switchBox == 'onn'){
    library(switchBox)
    P_dict <- SWAP.Filter.Wilcoxon(as.matrix(t(node.Data)), phenoGroup = as.factor(node.ClassLabels), featureNo = f_n)
    genePairs <- FindPairs(P_dict, node.Data, node.ClassLabels)
    print(P_dict)
    print(genePairs)
  }else{
  #First select the highly variable genes that correlate with the PCs
  P_dict <- FeatureSelector(Data = node.Data,
                             ClassLabels = node.ClassLabels,
                             num = 2000,
                             ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  #Then, select the genes as predictors if statistically DE between the classes.
  P_dict <- FeatureSelector2(Data = node.Data,
                            ClassLabels = node.ClassLabels,
                            num = f_n,
                            ...)
  }
  genePairs <- FindPairs(P_dict, node.Data, node.ClassLabels)

  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)
  if(labs_l[node] != "TaxaRoot"){
    childNodes <- GetChildNodeLeafs(tree = tree, node = node)
    childLeafs <- NULL
    for (i in which(lapply(childNodes, length) > 0)){
      childLeafs <- c(childLeafs, childNodes[i][[1]])
    }
    outGroupLeafs <- tree$tip.label[!tree$tip.label %in% childLeafs]
    node.outData <- droplevels(Rdata[which(Rdata$ClassLabels %in% outGroupLeafs), ])
    if(dim(node.outData)[1] > 500){
      node.outData <- node.outData[sample(rownames(node.outData), size = 500), ]#500 can be replaced with a better #
    }
    node.outData <- droplevels(subset(node.outData, select=c(P_dict)))
    node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
    node.Data <- rbind(node.Data, node.outData)
  }



  if(switchBox == 'on'){
    node.ClassLabels <- node.Data$ClassLabels
    node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))

    bins <- Binarize(genePairs, node.Data)
    node.Data <- as.data.frame(t(bins))
    node.Data$ClassLabels <- node.ClassLabels

  }
  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = TRUE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)
  node.mod <- caret::train(ClassLabels~.,
                           data = node.Data,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = FALSE,
                           proximity = FALSE,
                           outscale = FALSE,
                           preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)

  return(node.mod)
}

#' @param
#' @usage  genePairs <- FindPairs(P_dict, node.Data, node.ClassLabels)
FindPairs <- function(P_dict, node.Data, node.ClassLabels){
  library(tspair)
  genes <- P_dict
  df <- t(node.Data)
  #create the gene pairs:
  data <- df[P_dict, ]
  genePairs <- c()
  while(length(genes) >= 2){
    tsp1 <- tspcalc(as.matrix(data), node.ClassLabels)
    pairs <- rownames(tsp1$tspdat)
    if(!pairs[1] == pairs[2]){
      genePairs <- append(x = genePairs, paste(pairs[1], pairs[2], sep = "_"))
    }
    genes <- genes[!genes %in% pairs]
    data <- data[genes,]
    #print(paste(paste(pairs[1], pairs[2], sep = "_"), tsp1$tspscore, sep = " "))
  }
  return(genePairs)
}

#' @param
#' @usage  genePairs <- FindPairs(P_dict, node.Data, node.ClassLabels)
PairsPool <- function(P_dict, node.Data, node.ClassLabels, score_threshold=0.4){
  library(tspair)
  genes <- P_dict
  df <- t(node.Data)

  genePairs <- c()
  pairs.pool <- c()

  for(c in unique(as.character(node.ClassLabels))){
    genes <- P_dict
    df <- t(node.Data)
    #create the gene pairs:
    data <- df[P_dict, ]

    node.labs <- as.character(node.ClassLabels)
    node.labs[which(node.labs != c)] <- "A"
    node.labs[which(node.labs == c)] <- "B"

    while(length(genes) >= 2){
      tsp1 <- tspcalc(as.matrix(data), node.labs)
      pairs <- rownames(tsp1$tspdat)
      if((!pairs[1] == pairs[2]) && (tsp1$tspscore > score_threshold)){
        pairs.pool <- append(pairs.pool,pairs[1])
        pairs.pool <- append(pairs.pool,pairs[2])
      }
      genes <- genes[!genes %in% pairs]
      data <- data[genes,]
    }
  }
  return(unique(pairs.pool))
}


Binarize <- function(genePairs, node.Data){
  #Convert the node data to binary:
  df <- t(node.Data)
  bins <- matrix(0, nrow=length(genePairs), ncol=ncol(df))
  gPairs <- strsplit(genePairs, "_")
  for(i in seq(length(genePairs))){
    bins[i, ] <- as.numeric(df[gPairs[[i]][1], ] > df[gPairs[[i]][2], ])
  }
  rownames(bins) <- genePairs
  colnames(bins) <- colnames(df)
  return(bins)
}

RunGSEAforClusters <- function(SeuratObj, Cond1, Cond2, GeneSet=here('data/GeneSetDatabases/MousePath_GO_gmt.gmt'), outputDir=getwd(), ...) {

  SeuratObj <- SetAllIdent(SeuratObj, id = 'cluster.states')
  clusters <- sort(unique(as.numeric(SeuratObj@meta.data$tree.ident)))
  pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)

  for (i in clusters) {
    #Create GSEA input files
    filePrefix <- paste(Cond1,Cond2,"Cluster",i,"ExpMatrix",sep="_")
    CreateGSEAinput(SeuratObj = SeuratObj, Cond1 = Cond1, Cond2 = Cond2, clusterid = i, filePrefix = filePrefix, ...)
    #Run GSEA for the cluster i comparing cells from Cond1 and Cond2
    RunGSEA(InputPrefix = filePrefix, GeneSet=GeneSet, ...)

    setTxtProgressBar(pb,i)
    print(c('GSEA with Cluster', i,' is completed.'))
  }
}

RunDGEA <- function(SeuratObj, Cond1, Cond2, outputDir=getwd(), ...) {
  #First install MAST:
  #install.packages("BiocManager")
  #BiocManager::install("MAST")
  #library(scater)
  library(dplyr)
  SeuratObj <- SetAllIdent(SeuratObj, id = 'cluster.states')
  DEGs <- list()
  clusters <- sort(unique(as.numeric(SeuratObj@meta.data$tree.ident)))
  pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
  for (i in clusters) {
    ident.1 <- paste(Cond1,"-",i, sep = '')
    ident.2 <- paste(Cond2,'-',i, sep = '')
    a <- FindMarkers(SeuratObj, ident.1 = ident.1, ident.2 = ident.2, test.use = 'MAST', only.pos = FALSE, latent.vars = c('nUMI', 'percent.mito'), logfc.threshold = 0.1, min.pct = 0.05)
    a$gene <- rownames(a)
    a$cluster <- paste(i,sep = '')
    DEGs <- bind_rows(DEGs,a)
    setTxtProgressBar(pb,i)
    print(c('Finished with Cluster', i))
  }
  save(DEGs, file = paste(outputDir,"/",Cond1,"-",Cond2,"DEGs.Rdata",sep = ""))
  #unloadNamespace(ns = 'scater')
  return(DEGs)
}
#data_prep
CreateGSEAinput <- function(SeuratObj, Cond1, Cond2, outputDir=getwd(), clusterid, filePrefix, ...) {

  ident.1 <- paste(Cond1,"-",clusterid, sep = '')
  ident.2 <- paste(Cond2,'-',clusterid, sep = '')
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
  header1 <- "#1.2"
  header2 <- paste(dimensions[1],dimensions[2],sep="\t")
  write.table(header1, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE,col.names=FALSE,row.names=FALSE   )
  write.table(header2, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE  )
  write.table(new.df, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE, append=TRUE,col.names=NA   )

  #create a CLS file: for details: https://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS
  conditions <- c(ident.1,ident.2)
  header <- paste(dimensions[2], "2", "1",sep=" ")
  line2 <- paste("#",conditions[1],conditions[2], sep=" ")
  line3 <- paste( rep(c(conditions[1],conditions[2]), c(c1,c2)),sep = " " )
  write.table(header, file=paste(filePrefix,".cls",sep=""), sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
  write.table(line2,file=paste(filePrefix,".cls",sep=""), sep=" ", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE)
  linex <- line3[1]
  for (i in 2:length(line3)) {
    linex <- paste(linex,line3[i],sep =" ")}
  write.table(linex,file=paste(filePrefix,".cls",sep=""), sep=" ", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE)
}

RunGSEA <- function(InputPrefix, GeneSet, outputDir=getwd(), ...) {
  GSEA.program.location <- paste(srcdir,"/GSEA.1.1.R",sep="")
  source(GSEA.program.location, verbose=F, max.deparse.length=9999)

  doc.STRING <- paste(InputPrefix,substr(GeneSet,1,nchar(GeneSet)-4), sep="_")
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

SummarizeGSEAoutputs <- function(GSEAoutputDir="./") {
  #This function returns the table of all Enrichment results with corrected p-values.
  library(tidyverse)
  setwd(GSEAoutputDir)

  majorSummaryTable <- NULL
  GSreportsTable <- NULL
  mySumFiles <- list.files(pattern="*SUMMARY.RESULTS.REPORT*")

  for (i in 1:length(mySumFiles)) {

    sumTable <- read.delim(mySumFiles[i]) %>% as.tibble() %>% add_column(Comparison=strsplit(mySumFiles[i],"_Clust")[[1]][1],EnrichmentDirection_ClusterID=strsplit(mySumFiles[i],"\\.")[[1]][5])
    majorSummaryTable <- bind_rows(majorSummaryTable, sumTable)

    #for each Gene set, j, in the summary table:
    for(j in 1:length(read.delim(mySumFiles[i])[, 1])) {
      #the Gene Set j from the directory: Get the file prefix from the Summary file name + combine with gene set name + add ".txt" to the end.
      geneSetReportfile=list.files(pattern=paste(strsplit(mySumFiles[i], "\\.")[[1]][1], (read.delim(mySumFiles[i]) %>% as.tibble() %>% select(GS) %>% c())[[1]][j], "report", strsplit(mySumFiles[i], "\\.")[[1]][5], "*.txt", sep = "."))

      #if (!identical(geneSetReportfile, character(0))) {
      if (!identical(geneSetReportfile, character(0)) && (geneSetReportfile != "Subordinate_Control_Cluster_19_ExpMatrix_Calvin_manual_genesets.neuromuscular junction.report.Control-19.12.txt")) {

        gs.reporttable <-  read.delim(geneSetReportfile) %>%
          as.tibble() %>%
          dplyr::filter(CORE_ENRICHMENT == "YES") %>% # filter out genes which are not in the Leading Edge.
          add_column(
              Comparison = strsplit(mySumFiles[i], "_Clust")[[1]][1], #Create a column for Comparison type, ex; 'Dominant_Control'
              EnrichmentDirection_ClusterID = strsplit(mySumFiles[i], "\\.")[[1]][5], #Create a column for Enrichment direction, ex; 'Control-1'. This also shows the cluster id.
              GS = (read.delim(mySumFiles[i]) %>% as.tibble() %>% select(GS) %>% c())[[1]][j] #Create a column for Gene Set name.
            )
        GSreportsTable <- bind_rows(GSreportsTable, gs.reporttable)

      }else{break}#closes ifelse for report file existance.
    }#closes loop for j
  }#closes loop for i

  majorSummaryTable <- majorSummaryTable %>% as.tibble() %>% mutate(pAdj.Nom=p.adjust(NOM.p.val, method="BH")) %>% arrange(pAdj.Nom)
  sigtable <- majorSummaryTable %>% dplyr::filter(pAdj.Nom < 0.05) %>% unite(plotfilename, Comparison, GS, EnrichmentDirection_ClusterID, sep="*", remove = FALSE)
  #Write the main table and only significant enrichments to separate files:
  majorSummaryTable %>% write_tsv(file.path("All.Enrichment.stats.txt"))
  sigtable %>% write_tsv(file.path("significant.Enrichments.txt"))

  sig.GSreportsTable=NULL
  for(g in 1:dim(sigtable)[1]) {
    sig.g <- GSreportsTable %>% filter(Comparison == sigtable[g, ]$Comparison, EnrichmentDirection_ClusterID == sigtable[g, ]$EnrichmentDirection_ClusterID, GS == sigtable[g, ]$GS) %>% select(-SYMBOL, -DESC) %>% separate(EnrichmentDirection_ClusterID, into=c("EnrichmentDirection", "cluster"))
    sig.GSreportsTable <- bind_rows(sig.GSreportsTable, sig.g)
  }

  #return(majorSummaryTable)
  return(sig.GSreportsTable)
}

#' Get the score of the best class
#' @param vec vector of scores.
#' @param alphaList list of alpha threshold of each class.
GetScr <- function(vec, alphaList, ...){
  alphaList <- t(as.data.frame(alphaList))
  cert.vec <- vec > alphaList[, "Low.ext"]
  candits <- names(cert.vec)[cert.vec]
  if(identical(candits, character(0)) ){
    cls.max <- 0
  }else{
    cls.max <- max(vec[candits])
  }
  return(cls.max)
}


FeaturePrep <- function(SeuratObj, gene.set, ScoreName) {
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
#data_prep
SelectGenesBestLoadings <- function(trainingExpData, pcs, run.name="gene.select", num, caption="Highest loading") {
  #Second version:
  #p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom
  #Weighted gene picking depending on PC number: Initial PCs give more genes.
  #For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
  print("Performing PCA...")
  #Do PCA here on input data
  library(matrixStats)
  trainingExpData <- trainingExpData[, apply(trainingExpData, 2, var) != 0]
  trainingExpData <- trainingExpData[, which(colVars(as.matrix(trainingExpData)) > 0.05)]

  if (file.exists(paste(run.name,".train.prcomp.Rdata",sep=""))) {
    pcatrain <- get(load(paste(run.name,".train.prcomp.Rdata",sep="")))
  }else{
    pcatrain <- prcomp(trainingExpData, center = TRUE, scale=TRUE, rank. = pcs)
    save(pcatrain,file=paste(run.name,".train.prcomp.Rdata",sep=""))
  }

  print("Selecting the genes as best features...")

  trainingData <- get(load(paste(run.name,".trainingData.tmp.Rdata",sep = "")))
  pcadata <- data.frame(pcatrain$x, CellType = trainingData$CellType)
  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  pdf(paste(run.name,"_PCAplots.pdf",sep=""),width = 20,height = 10)
  pb <- txtProgressBar(min = 0, max = pcs, style = 3)
  for(i in 1:(pcs-1)) {
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[,i]),decreasing = TRUE),i]
    Gn <- round((pcs-i+1)*(num*2)/(pcs*(pcs+1)))
    TotalGenes <- as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai,Gn)),bestgenes=head(orderedpcai,Gn))
    load <- rbind(load, top)
    setTxtProgressBar(pb,i)
    pci <- paste("PC",i,sep="")
    pcj <- paste("PC",i+1,sep="")
    print(ggplot(pcadata, aes_string(x=pci, y=pcj, color="CellType"))+geom_point() )
    par()
    cat("\n",'Picking the best genes from first', i,'PCs is completed.',"\n")
  }
  dev.off()

  bestgenes <- unique(load$genes)
  return(bestgenes)
}
#data_prep
prepareDataset <- function(ExpressionData, CellLabels, run.name, PCs, featureGeneSet) {
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
  rm(trainingData)
  gc()

  if (missing(featureGeneSet)) {
    #Perform PCA on data (optional: only on training portion) and select genes for features
    genes <- as.character(SelectGenesBestLoadings(trainingExpData = train, run.name = run.name, pcs = PCs, num = 2000))
  }else{
    genes <- make.names(featureGeneSet)
  }#closes.if.missing.featureset

  trainingData <- get(load(paste(run.name,".trainingData.tmp.Rdata",sep = "")))

  trainingData.postPCA <- droplevels(trainingData[,c(genes,"CellType")])

  save(trainingData.postPCA, file=paste(run.name, ".trainingData.postPCA.data", sep = ""))
  file.remove(paste(run.name, ".trainingData.tmp.Rdata", sep = ""))

  return(trainingData.postPCA)
}

CellTyperTrainer2 <- function(ExpressionData, CellLabels, model.method="rf", run.name, do.splitTest=F, PCs, improve.rf=F, cv.k=5) {
  library(randomForest)
  #library(rfUtilities)
  library(tidyverse)
  library(caret)

  if (missing(PCs)) {
    PCs <- length(unique(CellLabels))
  }else{
    PCs <- PCs
  }

  ExpressionData <- as.matrix(ExpressionData)

  if (file.exists(paste(run.name,".trainingData.postPCA.data",sep = ""))) {
    print("Training data already exists...")
    trainingData <- get(load(paste(run.name,".trainingData.postPCA.data",sep = "")))
  }else{
    print("creating the training data...")
    trainingData <- prepareDataset(ExpressionData = ExpressionData, CellLabels = CellLabels, PCs = PCs, run.name = run.name)
  }

  #k-fold Cross Validation
  train.control <- trainControl(method="cv", number=cv.k, savePredictions = TRUE)
  model <- train(CellType~., data=trainingData, trControl=train.control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)
  save(model, file=paste(run.name,".RF_model_notImproved.Robj", sep = ""))

  if ((improve.rf == T) & (model.method == "rf")) {

    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))

    rfvotes <- as.data.frame(model$finalModel$votes)
    rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
    rfvotes$label <- apply(model$finalModel$votes, 1, function(x)  names(x)[which.max(x)] )
    rfvotes$inputnames <- rownames(rfvotes)
    e <- as.data.frame(table(rfvotes$label))
    Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()

    #Defaults
    filter.p <- 0.05
    bestvote.cutoff <- 0.7
    badinput.ids <- NULL
    round.n <- 1
    badinput.stats <- data.frame()
    currentscore <- Bestscore$median
    toss.n <- dim(rfvotes)[1]

    #While Median Best votes score overall is larger in the new iteration AND the number of inputs need to be tossed is larger %1 of all input data, then continue.
    #while (Bestscore$median >= currentscore && toss.n > round(0.01*dim(rfvotes)[1]) ) {
    while (round.n < 2 ) {#run this only once...

      print(paste("Round number ",round.n))
      print(paste("Current score is", currentscore,". toss.n is", toss.n, ". Fractions is", round(0.01*dim(rfvotes)[1])))

      for(i in 1:length(e$Var1)) {
        badinputs <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote.cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter.p*e$Freq[i]), wt=-bestvote) %>% select(inputnames) %>% c()
        badinputs.m <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote.cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter.p*e$Freq[i]), wt=-bestvote) %>% summarize(mean=mean(bestvote)) %>% c()
        badinput.ids <- unique(c(badinput.ids, badinputs$inputnames))
        classBestscore <- rfvotes %>% as.tibble() %>%  filter(., label == e$Var1[i]) %>% summarize(median=median(bestvote)) %>% c()
        class.badinput.stats <- data.frame(class=e$Var1[i], classInputSize=e$Freq[i], allBestscoreMedian=Bestscore$median, classBestscoreMedian=classBestscore$median, tossedInput=length(badinputs$inputnames), tossedBestvoteMean=badinputs.m$mean, iteration=round.n)
        badinput.stats <- rbind(badinput.stats, class.badinput.stats)
      }

      badinput.stats[is.nan(badinput.stats)] <- 0
      toss.n <- badinput.stats %>% as.tibble() %>% filter(., iteration == round.n) %>% summarise(n=sum(tossedInput)) %>% c()
      toss.n <- toss.n$n

      print(badinput.stats)

      #filter input using the bad input list generated in the previous iteration
      trainingData <- trainingData[which(!rownames(trainingData) %in% badinput.ids ),]

      #run the RF again with the updated training set
      model <- train(CellType~., data=trainingData, trControl=train.control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)

      #Check again to see if there is room to improve:
      rfvotes <- as.data.frame(model$finalModel$votes)
      rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
      rfvotes$label <- apply(model$finalModel$votes, 1, function(x)  names(x)[which.max(x)] )
      rfvotes$inputnames <- rownames(rfvotes)
      e <- as.data.frame(table(rfvotes$label))
      Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
      #update the round
      round.n = round.n + 1
    }#closes the while loop

  }#closes the improve option

  save(model, file=paste(run.name,".RF_model.Robj",sep = ""))
  return(model)
}

PlotPredictions2 <- function(SeuratObject, model, priorLabels, outputFilename="plotpredictions") {
  #Evaluate model prediction accuracy:
  conf.mat <- model$finalModel$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>% mutate(freq = 100*value/sum(value))
  class.n = length(model$finalModel$classes)
  errorSize <- as.data.frame(cbind(model$finalModel$confusion[,"class.error"], head(colSums(model$finalModel$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- getTrainPerf(model)["TrainAccuracy"]*100
  di <- round(sqrt(class.n),digits = 0)+1

  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +
    geom_tile(color = "white")+
    coord_equal()+
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

CellTyper2 <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions") {

  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
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

  if (!is.nan(mmDGm)) {
    if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
      warning("A significant portion of features are missing...")
      }
  }


  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class.n <- length(model$finalModel$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="raw")
  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

  if (missing(priorLabels)) {
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
      scale_x_discrete(limits = c("Prior", "Clusters", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
  }

  if (!missing(SeuratObject)) {
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

VarAsString <- function(x) {deparse(substitute(x))}

HTyper2 <- function(SeuratObject, testExpSet, models, priorLabels, outputFilename="plotpredictions") {

  #models is a list of of rf models
  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)

  colnames(testExpSet) <- make.names(colnames(testExpSet))

  Htable <- data.frame(cells = rownames(testExpSet))

  for (model in models) {
    modelname <- VarAsString(model)
    print(paste("Predicting with model",modelname,"...",sep = " "))

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

    #if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
    #  warning("A significant portion of features are missing...")
    #}

    rm(testsub, missingGenes, missingGenes.df)
    gc()
    #Predict
    library(entropy)
    testPred <- as.data.frame(predict(model, TestData, type = "prob"))
    class.n <- length(model$finalModel$classes)
    testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
    testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
    testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
    testPred$Prediction <- predict(model, TestData, type="raw")
    #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
    testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
    testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
    testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
    #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

    Htable <- data.frame(Htable, modelname = testPred$Prediction)

  }#closes models for loop

    if (missing(priorLabels)) {
      print("Prior class labels are not provided!")

    }else{
      #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
      priorLabels <- as.data.frame(priorLabels)
      colnames(priorLabels) <- c("Prior")
      Htable <- cbind(priorLabels,Htable)
      print(head(Htable))
      #Plot the crosscheck here:
      #Crosscheck Predictions
      library(tidyverse)
      library(alluvial)
      library(ggalluvial)
      Htable <- Htable %>% select(-cells) %>% droplevels()
      print(head(Htable))
      print(names(Htable))
      crx <- Htable %>% group_by_at(vars(one_of(names(Htable)))) %>% tally() %>% as.data.frame()
      print(head(crx))

      p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2], axis3 = names(crx)[3], axis4 = names(crx)[4] )) +
        geom_alluvium(aes_string(fill = names(crx)[4]), width = 1/12) +
        geom_stratum(width = 1/12, fill = "black", color = "red") +
        geom_label(stat = "stratum", label.strata = TRUE) +
        scale_x_discrete(limits = names(crx), expand = c(.05, .05)) +
        ggtitle("Predictions Cross-Check")

      save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
    }#closes missing PriorLabels

  if (!missing(SeuratObject)) {
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

PlotPredictions22 <- function(SeuratObject, outputFilename="plotpredictions") {

  #Prediction Performance
  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=FinalPrediction,fill=FinalPrediction),stat = "count")+
    geom_violin(aes(x=FinalPrediction,y=FinalBestProb*max(table(SeuratObject@meta.data$FinalPrediction)),fill=FinalPrediction))+
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$FinalPrediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))

  p4 <- TSNEPlot(SeuratObject, group.by="FinalPrediction",do.label=T, do.return = T)

  p3p4 <- cowplot::plot_grid(p3,p4, ncol = 1, nrow=2)
  save_plot(filename = paste(outputFilename,".prediction-stats.pdf",sep=""),plot = p3p4,base_height = 20, base_width = 24)

}

HTyper22 <- function(SeuratObject, testExpSet, taxTable, models, priorLabels, outputFilename="plotpredictions") {

  #models is a list of of rf models
  library(caret)
  library(randomForest)
  library(tidyverse)
  if (missing(taxTable)) {stop("Please provide a proper taxanomy table as a dataframe with 'taxTable' ... exiting!")}else{taxtable <- taxTable}

  if (!missing(SeuratObject)) {
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
    }#Closes missing(SeuratObj)

  colnames(testExpSet) <- make.names(colnames(testExpSet))

  Htable <- data.frame(cells = rownames(testExpSet))
  Bigtable <- data.frame(cells = rownames(testExpSet))
  i <- 1
  for (model in models) {
    modelname <- VarAsString(model)
    print(paste("Predicting with model",modelname,"...",sep = " "))

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

    #if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
    #  warning("A significant portion of features are missing...")
    #}

    rm(testsub, missingGenes, missingGenes.df)
    gc()
    #Predict
    library(entropy)
    testPred <- as.data.frame(predict(model, TestData, type = "prob"))
    class.n <- length(model$finalModel$classes)
    testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
    testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
    testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
    testPred$Prediction <- predict(model, TestData, type="raw")
    #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
    testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
    testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
    testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
    #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

    colnames(testPred) <- paste(modelname,i, names(testPred),sep = ".")
    Bigtable <- cbind(Bigtable, testPred)
    Htable <- data.frame(Htable, modelname = testPred[,  paste(modelname,i,"Prediction",sep = ".")])

    i=i+1
  }#closes models for loop


  ConditionalProbTable <- data.frame(matrix(ncol = 0, nrow = length(rownames(Bigtable))))
  leafNames <- NULL
  #For each leaf:
  for(j in 1:dim(taxtable)[1]) {
    leafName <- paste(taxtable[j,dim(taxtable)[2]],sep = "")
    print(leafName)
    #Calculate the Conditional Probabilities
    nodeNames <- NULL
    for(i in 1:length(models.list)) {
      print(paste("model",i, taxtable[j,i], sep="."))
      nodeNames <- c(nodeNames, paste("model",i, taxtable[j,i], sep="."))
    }
    ConditionalProbTable <- cbind(ConditionalProbTable, matrixStats::rowProds(as.matrix(Bigtable[,nodeNames])))
    leafNames <- c(leafNames, leafName)
  }
  colnames(ConditionalProbTable) <- leafNames
  ConditionalProbTable$FinalBestProb <- apply(ConditionalProbTable, 1, function(x) max(x) )
  ConditionalProbTable$FinalPrediction <- colnames(ConditionalProbTable[,which(!colnames(ConditionalProbTable) %in% c("FinalBestProb"))])[apply(ConditionalProbTable[,which(!colnames(ConditionalProbTable) %in% c("FinalBestProb"))],1,which.max)]

  Htable$FinalPrediction <- ConditionalProbTable$FinalPrediction


  if (missing(priorLabels)) {print("Prior class labels are not provided!")}else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    Htable <- cbind(priorLabels,Htable)
    print(head(Htable))
    save(Htable, file=paste(outputFilename,".prediction-crosscheck.Full_htable.Rdata",sep=""))

    #Plot the crosscheck here: For interactive Sankey diagram: https://www.r-graph-gallery.com/sankey-diagram/
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    Htable <- Htable %>% select(-cells) %>% droplevels()
    print(head(Htable))
    print(names(Htable))
    crx <- Htable %>% group_by_at(vars(one_of(names(Htable)))) %>% tally() %>% as.data.frame()
    print(head(crx))
    crx.f <- crx %>% mutate(freq = n*100 / sum(n)) %>% filter(freq > 1)

    p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2], axis3 = names(crx)[3], axis4 = names(crx)[4], axis5 = names(crx)[5], axis6 = names(crx)[6] )) +
      geom_alluvium(aes_string(fill = names(crx)[6]), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("PriorLabels","Zeisel.Tax.Rank1","Zeisel.Tax.Rank2","Zeisel.Tax.Rank3","Zeisel.Tax.Rank4","FinalPrediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    p6 <- ggplot(crx.f,aes_string(y = "n", axis1 = names(crx.f)[1], axis2 = names(crx.f)[2], axis3 = names(crx.f)[3], axis4 = names(crx.f)[4], axis5 = names(crx.f)[5], axis6 = names(crx.f)[6])) +
      geom_alluvium(aes_string(fill = names(crx.f)[6]), width = 0, knot.pos = 1/4) +
      guides(fill = FALSE)+
      geom_stratum(width = 1/12, fill = "grey", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE ) +
      ylab("Frequency")+
      scale_x_discrete(limits = c("PriorLabels","Zeisel.Tax.Rank1","Zeisel.Tax.Rank2","Zeisel.Tax.Rank3","Zeisel.Tax.Rank4","FinalPrediction"), expand = c(.05, .05)) +
      ggtitle(paste("Predictions","Cross-Check",sep = " "))

    pdf(paste(outputFilename,".prediction-crosscheck.pdf",sep=""),width = 40,height = 30)
    print(p5)
    print(p6)
    dev.off()
    save(crx, file=paste(outputFilename,".prediction-crosscheck.htable.Rdata",sep=""))
  }#closes missing PriorLabels

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(ConditionalProbTable))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, ConditionalProbTable)

    PlotPredictions22(SeuratObject = SeuratObject, outputFilename = outputFilename)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(Bigtable) <- rownames(testExpSet)
    return(Bigtable)
  }#Closes missing(SeuratObj)
}#closes the function

TrainPrep <- function(model, modeSeurat, RankLabellist ) {
  imp <- as.data.frame(model$finalModel$importance)
  features <- rownames(head(imp[order(imp$MeanDecreaseGini, decreasing=T),],200))
  cells <- rownames(model$finalModel$votes)
  RLR <- modeSeurat@meta.data[,RankLabellist]
  colns <- NULL
  for(i in 1:(length(RankLabellist)-1 ) ) {colns <- c(colns,c(paste("R",i,sep="")))}
  print(colns)
  colnames(RLR) <- c(colns,"CellType")
  trainingData <- as.data.frame(t(as.matrix(modeSeurat@data)))
  #It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
  names(trainingData) <- make.names(names(trainingData))
  trainingData <- cbind(trainingData[,features], RLR)
  #Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
  trainingData$CellType <- factor(trainingData$CellType)
  gc()
  return(trainingData)
}

RecursiveTrainer <- function(trainingData, model.method="rf", run.name, cv.k=5) {
  library(randomForest)
  library(rfUtilities)
  library(tidyverse)
  library(caret)

  #k-fold Cross Validation
  #train.control <- trainControl(method="cv", number=cv.k, savePredictions = TRUE)
  model <- train(CellType~., data=trainingData, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)

  save(model, file=paste(run.name,".RFrec_model.Robj", sep = ""))

  return(model)
}

####################################Previous version functions###############################################
CellTyperTrainer <- function(ExpressionData, CellLabels, run.name, do.splitTest=F, PCs, improve=T) {
  library(randomForest)
  library(rfUtilities)
  library(tidyverse)

  if (missing(PCs)) {
    PCs=length(unique(CellLabels))
  }else{
    PCs=PCs
  }

  ExpressionData <- as.matrix(ExpressionData)

  if (file.exists(paste(run.name,".trainingData.postPCA.data",sep = ""))) {
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

  if (improve == T) {
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))

    rfvotes <- as.data.frame(rf$votes)
    rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
    rfvotes$label <- apply(rf$votes, 1, function(x)  names(x)[which.max(x)] )
    rfvotes$inputnames <- rownames(rfvotes)
    e <- as.data.frame(table(rfvotes$label))
    Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()

    #Defaults
    filter.p = 0.05
    bestvote.cutoff = 0.7
    badinput.ids <- NULL
    round.n=1
    badinput.stats <- data.frame()
    currentscore = Bestscore$median
    toss.n = dim(rfvotes)[1]

    #While Median Best votes score overall is larger in the new iteration AND the number of inputs need to be tossed is larger %1 of all input data, then continue.
    #while (Bestscore$median >= currentscore && toss.n > round(0.01*dim(rfvotes)[1]) ) {
    while (round.n < 2 ) {#run this only once...

      print(paste("Round number ",round.n))
      print(paste("Current score is", currentscore,". toss.n is", toss.n, ". Fractions is", round(0.01*dim(rfvotes)[1])))

      for(i in 1:length(e$Var1)) {
        badinputs <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote.cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter.p*e$Freq[i]), wt=-bestvote) %>% select(inputnames) %>% c()
        badinputs.m <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote.cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter.p*e$Freq[i]), wt=-bestvote) %>% summarize(mean=mean(bestvote)) %>% c()
        badinput.ids <- unique(c(badinput.ids, badinputs$inputnames))
        classBestscore <- rfvotes %>% as.tibble() %>%  filter(., label == e$Var1[i]) %>% summarize(median=median(bestvote)) %>% c()
        #print(paste("The number of input dropped for",e$Var1[i],"class is", length(badinputs$inputnames),sep = " "))
        class.badinput.stats <- data.frame(class=e$Var1[i], classInputSize=e$Freq[i], allBestscoreMedian=Bestscore$median, classBestscoreMedian=classBestscore$median, tossedInput=length(badinputs$inputnames), tossedBestvoteMean=badinputs.m$mean, iteration=round.n)
        badinput.stats <- rbind(badinput.stats, class.badinput.stats)
      }

      badinput.stats[is.nan(badinput.stats)] <- 0
      toss.n <- badinput.stats %>% as.tibble() %>% filter(., iteration == round.n) %>% summarise(n=sum(tossedInput)) %>% c()
      toss.n <- toss.n$n

      print(badinput.stats)

      #filter input using the bad input list generated in the previous iteration
      trainingData <- trainingData[which(!rownames(trainingData) %in% badinput.ids ),]

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
      round.n = round.n + 1
    }#closes the while loop

  }#closes the improve option

  save(rf, file=paste(run.name,".RF_model.Robj",sep = ""))
  return(rf)
}

PlotPredictions <- function(SeuratObject, model, save.pdf=T, outputFilename="plotpredictions") {
  #Evaluate model prediction accuracy:
  conf.mat <- model$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>%
    mutate(freq = 100*value/sum(value))

  class.n = length(model$classes)

  pdf(paste(outputFilename,".pdf",sep=""),width= 1.5*class.n, height = 1.5*class.n)


  require(gridExtra)

  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +
    geom_tile(color = "white")+
    coord_equal()+
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

CellTyper <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions") {

  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
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
  class.n <- length(model$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="response")

  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

  if (missing(priorLabels)) {
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

  if (!missing(SeuratObject)) {

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

HTyper <- function(SeuratObject, testExpSet, models, priorLabels, outputFilename="plotpredictions") {

  #models is a list of of rf models
  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)

  colnames(testExpSet) <- make.names(colnames(testExpSet))

  Htable <- data.frame(cells = rownames(testExpSet))

  for (model in models) {
    modelname <- VarAsString(model)
    print(paste("Predicting with model",modelname,"...",sep = " "))

    #Prepare Test Expression set
    testsub <- testExpSet[,which(colnames(testExpSet) %in% attributes(model$terms)$term.labels)]
    missingGenes <- attributes(model$terms)$term.labels[which(!attributes(model$terms)$term.labels %in% colnames(testExpSet))]
    print(model$importance[missingGenes,])

    missingGenes.df <- data.frame(matrix(0, ncol = length(missingGenes), nrow = length(rownames(testExpSet))))
    colnames(missingGenes.df) <- missingGenes
    TestData <- cbind(testsub, missingGenes.df)
    TestData <- TestData[,attributes(model$terms)$term.labels]
    mmDGm <- mean(model$importance[missingGenes,])
    mmDGf <- mean(model$importance[which(!attributes(model$terms)$term.labels %in% missingGenes),])

    cat("Number of Features (genes) to be considered is", length(colnames(testsub)), '\n',
        "Number of missing Features set to zero is", length(missingGenes), '\n',
        "Mean MeanDecreaseGini of missing genes is", mmDGm, '\n',
        "Mean MeanDecreaseGini of Featured genes is", mmDGf, '\n',
        sep = ' ')

    #if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
    #  warning("A significant portion of features are missing...")
    #}

    rm(testsub, missingGenes, missingGenes.df)
    gc()
    #Predict
    library(entropy)
    testPred <- as.data.frame(predict(model, TestData, type = "prob"))
    class.n <- length(model$classes)
    testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
    testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
    testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
    testPred$Prediction <- predict(model, TestData, type="response")
    #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
    testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
    testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.05) | (Diff <= 0.05), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
    testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.05) | (Diff <= 0.05), "Undetermined", "Detected")) %>% as.data.frame()
    #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

    Htable <- data.frame(Htable, modelname = testPred$Prediction)

  }#closes models for loop

  if (missing(priorLabels)) {
    print("Prior class labels are not provided!")

  }else{
    # Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    Htable <- cbind(priorLabels,Htable)
    print(head(Htable))
    #Plot the crosscheck here:
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    Htable <- Htable %>% select(-cells) %>% droplevels()
    print(head(Htable))
    print(names(Htable))
    crx <- Htable %>% group_by_at(vars(one_of(names(Htable)))) %>% tally() %>% as.data.frame()
    print(head(crx))

    p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2], axis3 = names(crx)[3], axis4 = names(crx)[4], axis5 = names(crx)[5] )) +
      geom_alluvium(aes_string(fill = names(crx)[5]), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = names(crx), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
    save(crx, file=paste(outputFilename,".prediction-crosscheck.htable.Rdata",sep=""))
  }  # closes missing PriorLabels

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)

    #PlotPredictions(SeuratObject = SeuratObject, model = model, outputFilename = outputFilename)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }  # Closes missing(SeuratObj)
}  # closes the function
