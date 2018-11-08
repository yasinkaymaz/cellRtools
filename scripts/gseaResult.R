
GSEA.ConsPlot <- function(V, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {

# Plots a heatmap plot of a consensus matrix

     cols <- length(V[1,])
     B <- matrix(0, nrow=cols, ncol=cols)
     max.val <- max(V)
     min.val <- min(V)
     for (i in 1:cols) {
         for (j in 1:cols) {
             k <- cols - i + 1
	     B[k, j] <-  max.val - V[i, j] + min.val
          }
     }

 

#     col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
     col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))

#     max.size <- max(nchar(col.names))
     par(mar = c(5, 15, 15, 5))
     image(1:cols, 1:cols, t(B), col = col.map, axes=FALSE, main=main, sub=sub, xlab= xlab, ylab=ylab)

     for (i in 1:cols) {
         col.names[i]  <- substr(col.names[i], 1, 25)
     }
     col.names2 <- rev(col.names)

     size.col.char <- ifelse(cols < 15, 1, sqrt(15/cols))

     axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.col.char, font.axis=1, line=-1)
     axis(3, at=1:cols, labels=col.names, adj= 1, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=1, line=-1)

     return()
}

GSEA.HeatMapPlot2 <- function(V, row.names = "NA", col.names = "NA", main = " ", sub = " ", xlab=" ", ylab=" ", color.map = "default") {
#
# Plots a heatmap of a matrix

       n.rows <- length(V[,1])
       n.cols <- length(V[1,])

       if (color.map == "default") {
         color.map <- rev(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5))
       }

        heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
        heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]

        par(mar = c(7, 15, 5, 5))
        image(1:n.cols, 1:n.rows, t(heatm), col=color.map, axes=FALSE, main=main, sub = sub, xlab= xlab, ylab=ylab)

        if (length(row.names) > 1) {
            size.row.char <- ifelse(n.rows < 15, 1, sqrt(15/n.rows))
            size.col.char <- ifelse(n.cols < 15, 1, sqrt(10/n.cols))
#            size.col.char <- ifelse(n.cols < 2.5, 1, sqrt(2.5/n.cols))
            for (i in 1:n.rows) {
               row.names[i] <- substr(row.names[i], 1, 40)
            }
            row.names <- row.names[seq(n.rows, 1, -1)]
            axis(2, at=1:n.rows, labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=1, line=-1)
        }

        if (length(col.names) > 1) {
           axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
        }

	return()
}





GSEA.write.gct <- function (gct, filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct)[1], "\t", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")
    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")
    names <- names(gct)
    cat("\t", names[1], file = f, append = TRUE, sep = "")
    for (j in 2:length(names)) {
        cat("\t", names[j], file = f, append = TRUE, sep = "")
    }
    cat("\n", file = f, append = TRUE, sep = "\t")
    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] +  2)
    m[, 1] <- row.names(gct)
    m[, 2] <- row.names(gct)
    index <- 3
    for (i in 1:dim(gct)[2]) {
        m[, index] <- gct[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)
    return(gct)
}




###########################################
GSEA.Analyze.Sets <- function(
  directory,
  topgs = "",
  non.interactive.run = F,
  height = 12,
  width = 17) {
  
  file.list <- list.files(directory)
#print(directory)
  files <- file.list[regexpr(pattern = ".report.", file.list) > 1]
#print(files)
  max.sets <- length(files)
  
  set.table <- matrix(nrow = max.sets, ncol = 5)
  
  for (i in 1:max.sets) {
    temp1 <-  strsplit(files[i], split=".report.")
#print(temp1[[1]][1])
#    temp2 <-  strsplit(temp1[[1]][1], split=".")
    temp2 <-  strsplit(sub("\\.","*", temp1[[1]][1]),"\\*" )
#print(paste("temp2:",temp2))
    s <- length(temp2[[1]])
#print(temp2[[1]][2])
#    prefix.name <- paste(temp2[[1]][1:(s-1)], sep="", collapse="")
    prefix.name <- paste(temp2[[1]][1], sep="", collapse="")
#print(prefix.name)
    set.name <- temp2[[1]][s]
#    temp3 <-  strsplit(temp1[[1]][2], split=".")
    temp3 <-  strsplit(sub("\\.","*", temp1[[1]][2]),"\\*" )
#print(temp3[[1]][2])
    phenotype <- temp3[[1]][1]
#    seq.number <- temp3[[1]][2]
    seq.number <- strsplit(temp3[[1]][2], split=".txt")
#print(seq.number)
#    dataset <- paste(temp2[[1]][1:(s-1)], sep="", collapse=".")
    dataset <- paste(temp2[[1]][1], sep="", collapse="\\.")
#print(dataset)    
    set.table[i, 1] <- files[i]
    
    set.table[i, 3] <- phenotype
    set.table[i, 4] <- as.numeric(seq.number)
    set.table[i, 5] <- dataset
#print(set.table)    
    #      set.table[i, 2] <- paste(set.name, dataset, sep ="", collapse="")
    set.table[i, 2] <- substr(set.name, 1, 20) 
  }
#  print("Loop closed!")
#print(prefix.name)
  print(c("set name=", prefix.name))
  doc.string <- prefix.name
  
  set.table <- noquote(set.table)

  phen.order <- order(set.table[, 3], decreasing = T)
  set.table <- set.table[phen.order,]
#print(set.table)    
  phen1 <- names(table(set.table[,3]))[1]
  phen2 <- names(table(set.table[,3]))[2]
#print(phen2)
  set.table.phen1 <- set.table[set.table[,3] == phen1,]

  set.table.phen2 <- set.table[set.table[,3] == phen2,]
#print(set.table.phen2)
#print(set.table.phen1)  
  seq.order <- order(as.numeric(set.table.phen1[, 4]), decreasing = F)
  set.table.phen1 <- set.table.phen1[seq.order,]
  seq.order <- order(as.numeric(set.table.phen2[, 4]), decreasing = F)
  set.table.phen2 <- set.table.phen2[seq.order,]
  
  #   max.sets.phen1 <- length(set.table.phen1[,1])
  #   max.sets.phen2 <- length(set.table.phen2[,1])
  
  if (topgs == "") {
 
    max.sets.phen1 <- length(set.table.phen1[,1])
    max.sets.phen2 <- length(set.table.phen2[,1])
  } else {

    max.sets.phen1 <- ifelse(topgs > length(set.table.phen1[,1]), length(set.table.phen1[,1]), topgs) 
    max.sets.phen2 <- ifelse(topgs > length(set.table.phen2[,1]), length(set.table.phen2[,1]), topgs)
  }

  # Analysis for phen1
  
  leading.lists <- NULL
  for (i in 1:max.sets.phen1) {
    inputfile <- paste(directory, set.table.phen1[i, 1], sep="", collapse="")
#print(inputfile)
    gene.set <- read.table(file=inputfile, sep="\t", header=T, comment.char="", as.is=T)
    leading.set <- as.vector(gene.set[gene.set[,"CORE_ENRICHMENT"] == "YES", "SYMBOL"])
    leading.lists <- c(leading.lists, list(leading.set))
    if (i == 1) {
      all.leading.genes <- leading.set 
    } else{
      all.leading.genes <- union(all.leading.genes, leading.set)
    }
  }
  max.genes <- length(all.leading.genes)
  M <- matrix(0, nrow=max.sets.phen1, ncol=max.genes)
  for (i in 1:max.sets.phen1) {
    M[i,] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch=0))   # notice that the sign is 0 (no tag) or 1 (tag) 
  }
  
  Inter <- matrix(0, nrow=max.sets.phen1, ncol=max.sets.phen1)
  for (i in 1:max.sets.phen1) {
    for (j in 1:max.sets.phen1) {
      Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]]))/length(union(leading.lists[[i]], leading.lists[[j]]))
    }
  }
  
  Itable <- data.frame(Inter)
print(set.table.phen1)
#print(set.table.phen1[1:max.sets.phen1, 2])
#  names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
  names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
#  row.names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
   row.names(Itable) <- set.table.phen1[1:max.sets.phen1, 1] 
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, sep="", collapse="")
      windows(height = width, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = width, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  }
  
  GSEA.ConsPlot(Itable, col.names = set.table.phen1[1:max.sets.phen1, 2], main = " ", sub=paste("Leading Subsets Overlap ", doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ")
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  # Save leading subsets in a GCT file
  
  D.phen1 <- data.frame(M)
  names(D.phen1) <- all.leading.genes
#  row.names(D.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
  row.names(D.phen1) <- set.table.phen1[1:max.sets.phen1, 1]
  output <- paste(directory, doc.string, ".leading.genes.", phen1, ".gct", sep="")
  GSEA.write.gct(D.phen1, filename=output)
  
  # Save leading subsets as a single gene set in a .gmt file
  
  row.header <- paste(doc.string, ".all.leading.genes.", phen1, sep="")
  output.line <- paste(all.leading.genes, sep="\t", collapse="\t")
  output.line <- paste(row.header, row.header, output.line, sep="\t", collapse="")
  output <- paste(directory, doc.string, ".all.leading.genes.", phen1, ".gmt", sep="")
  write(noquote(output.line), file = output, ncolumns = length(output.line))
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, sep="", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  cmap <-  c("#AAAAFF", "#111166")
  GSEA.HeatMapPlot2(V = data.matrix(D.phen1), row.names = row.names(D.phen1), col.names = names(D.phen1), main = "Leading Subsets Assignment", sub = paste(doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  DT1.phen1 <- data.matrix(t(D.phen1))
  DT2.phen1 <- data.frame(DT1.phen1)
  names(DT2.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
  row.names(DT2.phen1) <- all.leading.genes
  #   GSEA.write.gct(DT2.phen1, filename=outputfile2.phen1)
  
  # Analysis for phen2
  
  leading.lists <- NULL
  for (i in 1:max.sets.phen2) {
    inputfile <- paste(directory, set.table.phen2[i, 1], sep="", collapse="")
    gene.set <- read.table(file=inputfile, sep="\t", header=T, comment.char="", as.is=T)
    leading.set <- as.vector(gene.set[gene.set[,"CORE_ENRICHMENT"] == "YES", "SYMBOL"])
    leading.lists <- c(leading.lists, list(leading.set))
    if (i == 1) {
      all.leading.genes <- leading.set 
    } else{
      all.leading.genes <- union(all.leading.genes, leading.set)
    }
  }
  max.genes <- length(all.leading.genes)
  M <- matrix(0, nrow=max.sets.phen2, ncol=max.genes)
  for (i in 1:max.sets.phen2) {
    M[i,] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch=0))   # notice that the sign is 0 (no tag) or 1 (tag) 
  }
  
  Inter <- matrix(0, nrow=max.sets.phen2, ncol=max.sets.phen2)
  for (i in 1:max.sets.phen2) {
    for (j in 1:max.sets.phen2) {
      Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]]))/length(union(leading.lists[[i]], leading.lists[[j]]))
    }
  }
  
  Itable <- data.frame(Inter)
  names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
#  row.names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(Itable) <- set.table.phen2[1:max.sets.phen2, 1]
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, sep="", collapse="")
      windows(height = width, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = width, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  }
  
  GSEA.ConsPlot(Itable, col.names = set.table.phen2[1:max.sets.phen2, 2], main = " ", sub=paste("Leading Subsets Overlap ", doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ")
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  # Save leading subsets in a GCT file
  
  D.phen2 <- data.frame(M)
  names(D.phen2) <- all.leading.genes
#  row.names(D.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(D.phen2) <- set.table.phen2[1:max.sets.phen2, 1]
  output <- paste(directory, doc.string, ".leading.genes.", phen2, ".gct", sep="")
  GSEA.write.gct(D.phen2, filename=output)
  
  # Save primary subsets as a single gene set in a .gmt file
  
  row.header <- paste(doc.string, ".all.leading.genes.", phen2, sep="")
  output.line <- paste(all.leading.genes, sep="\t", collapse="\t")
  output.line <- paste(row.header, row.header, output.line, sep="\t", collapse="")
  output <- paste(directory, doc.string, ".all.leading.genes.", phen2, ".gmt", sep="")
  write(noquote(output.line), file = output, ncolumns = length(output.line))
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, sep="", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  cmap <-  c("#AAAAFF", "#111166")
  GSEA.HeatMapPlot2(V = data.matrix(D.phen2), row.names = row.names(D.phen2), col.names = names(D.phen2), main = "Leading Subsets Assignment", sub = paste(doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  DT1.phen2 <- data.matrix(t(D.phen2))
  DT2.phen2 <- data.frame(DT1.phen2)
  names(DT2.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(DT2.phen2) <- all.leading.genes
  #   GSEA.write.gct(DT2.phen2, filename=outputfile2.phen2)
  
  # Resort columns and rows for phen1
  
  A <- data.matrix(D.phen1)
  A.row.names <- row.names(D.phen1)
  A.names <- names(D.phen1)
  
  # Max.genes
  
  #   init <- 1
  #   for (k in 1:max.sets.phen1) { 
  #      end <- which.max(cumsum(A[k,]))
  #      if (end - init > 1) {
  #         B <- A[,init:end]
  #         B.names <- A.names[init:end]
  #        dist.matrix <- dist(t(B))
  #         HC <- hclust(dist.matrix, method="average")
  ##         B <- B[,HC$order] + 0.2*(k %% 2)
  #        B <- B[,HC$order] 
  #         A[,init:end] <- B
  #         A.names[init:end] <- B.names[HC$order]
  #         init <- end + 1
  #     }
  #   }
  
  #   windows(width=14, height=10)
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = "  ", main = paste("Primary Sets Assignment - ", doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ") 
  
  dist.matrix <- dist(t(A))
  HC <- hclust(dist.matrix, method="average")
  A <- A[, HC$order]
  A.names <- A.names[HC$order]
  
  dist.matrix <- dist(A)
  HC <- hclust(dist.matrix, method="average")
  A <- A[HC$order,]
  A.row.names <- A.row.names[HC$order]
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, sep="", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  
  
  cmap <-  c("#AAAAFF", "#111166")
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  GSEA.HeatMapPlot2(V = t(A), row.names = A.names, col.names = A.row.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  text.filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".txt", sep="", collapse="")
  line.list <- c("Gene", A.row.names)
  line.header <- paste(line.list, collapse="\t")
  line.length <- length(A.row.names) + 1
  write(line.header, file = text.filename, ncolumns = line.length)
  write.table(t(A), file=text.filename, append = T, quote=F, col.names= F, row.names=T, sep = "\t")
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  
  
  
  
  
  # resort columns and rows for phen2
  
  A <- data.matrix(D.phen2)
  A.row.names <- row.names(D.phen2)
  A.names <- names(D.phen2)
  
  # Max.genes
  
  #   init <- 1
  #   for (k in 1:max.sets.phen2) { 
  #      end <- which.max(cumsum(A[k,]))
  #      if (end - init > 1) {
  #         B <- A[,init:end]
  #         B.names <- A.names[init:end]
  #         dist.matrix <- dist(t(B))
  #         HC <- hclust(dist.matrix, method="average")
  ##         B <- B[,HC$order] + 0.2*(k %% 2)
  #        B <- B[,HC$order] 
  #         A[,init:end] <- B
  #         A.names[init:end] <- B.names[HC$order]
  #         init <- end + 1
  #     }
  #   }
  
  #  windows(width=14, height=10)
  #  GESA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = "  ", main = paste("Primary Sets Assignment - ", doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ") 
  
  dist.matrix <- dist(t(A))
  HC <- hclust(dist.matrix, method="average")
  A <- A[, HC$order]
  A.names <- A.names[HC$order]
  
  dist.matrix <- dist(A)
  HC <- hclust(dist.matrix, method="average")
  A <- A[HC$order,]
  A.row.names <- A.row.names[HC$order]
  
  if (non.interactive.run == F) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, sep="", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep="", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  cmap <-  c("#AAAAFF", "#111166")
  
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  GSEA.HeatMapPlot2(V = t(A), row.names =A.names , col.names = A.row.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  text.filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".txt", sep="", collapse="")
  line.list <- c("Gene", A.row.names)
  line.header <- paste(line.list, collapse="\t")
  line.length <- length(A.row.names) + 1
  write(line.header, file = text.filename, ncolumns = line.length)
  write.table(t(A), file=text.filename, append = T, quote=F, col.names= F, row.names=T, sep = "\t")
  
  if (non.interactive.run == F) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
}



#GSEA.Analyze.Sets(
#  directory = "/home/yasinkaymaz/Downloads/top.tmp/",
#  topgs = 20,
#  height = 16,
#  width = 16
#)

