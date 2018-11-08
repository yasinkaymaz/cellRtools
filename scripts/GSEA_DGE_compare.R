#

load("~/Downloads/mn.RObj")
mn = SetAllIdent(mn, id = 'cluster.states')
clusters = sort(unique(as.numeric(mn@meta.data$tree.ident)))

DvsC.DEGs = list()
pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
for (i in clusters) {
  ident.1 = paste("Dominant-",i, sep = '')
  ident.2 = paste('Control-',i, sep = '')
  a =FindMarkers(mn, ident.1 = ident.1, ident.2 = ident.2, test.use = 'negbinom', 
                 only.pos = FALSE, min.cells.group = 1,
                 latent.vars = c('nUMI','percent.mito','batch'))
  a$gene = rownames(a)
  a$cluster = paste(i,sep = '')
  DvsC.DEGs = bind_rows(DvsC.DEGs,a)
  setTxtProgressBar(pb,i)
  print(c('Finished with Cluster', i))
}

# Subordinate vs control
SvsC.DEGs = list()
pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
for (i in clusters) {
  ident.1 = paste("Subordinate-",i, sep = '')
  ident.2 = paste('Control-',i, sep = '')
  a =FindMarkers(mn, ident.1 = ident.1, ident.2 = ident.2,  test.use = 'negbinom', 
                 only.pos = FALSE,  min.cells.group = 1,
                 latent.vars = c('nUMI','percent.mito','batch'))
  a$gene = rownames(a)
  a$cluster = paste(i,sep = '')
  SvsC.DEGs = bind_rows(SvsC.DEGs,a)
  setTxtProgressBar(pb,i)
  print(c('Finished with Cluster', i))
}

# Dominant vs Subordinate
DvsS.DEGs = list()
pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
for (i in clusters) {
  ident.1 = paste("Dominant-",i, sep = '')
  ident.2 = paste('Subordinate-',i, sep = '')
  a =FindMarkers(mn, ident.1 = ident.1, ident.2 = ident.2, test.use = 'negbinom', 
                 only.pos = FALSE, min.cells.group = 1,
                 latent.vars = c('nUMI','percent.mito','batch'))
  a$gene = rownames(a)
  a$cluster = paste(i,sep = '')
  DvsS.DEGs = bind_rows(DvsS.DEGs,a)
  setTxtProgressBar(pb,i)
  print(c('Finished with Cluster',i))
}

# combine into one dataframe and filter
DEGs = dplyr::bind_rows(DvsC.DEGs,DvsS.DEGs,SvsC.DEGs, .id = 'id')
DEGs$fdr = p.adjust(DEGs$p_val, method = 'fdr')
DEGs$id = factor(DEGs$id)
levels(DEGs$id)= list('DvsC'= '1' ,'DvsS' ='2','SvsC' = '3')
DEGs$cluster = as.numeric(DEGs$cluster)
save(DEGs, file='/Users/yasinkaymaz/Documents/DulacLab/Eric-Adam/gseaTest/MDT_mn/allDEGs.Rdata')
#load('/Users/yasinkaymaz/Documents/DulacLab/Eric-Adam/gseaTest/MDT_mn/allDEGs.Rdata')

sigEnGS <- read.delim("/Users/yasinkaymaz/Documents/DulacLab/Eric-Adam/gseaTest/MDT_mn/significant.Enrichments.txt",header=TRUE)
setwd("/Users/yasinkaymaz/Documents/DulacLab/Eric-Adam/gseaTest/MDT_mn/")
library(tibble)

paste(sigEnGS[1,]$Comparison,"_","Cluster_",str_split(sigEnGS[i,]$EnrichmentDirection_ClusterID,"-")[[1]][2],"_ExpMatrix_Calvin_manual_genesets.",sigEnGS[1,]$GS,".report.",sigEnGS[1,]$EnrichmentDirection_ClusterID,"*.txt",sep = "")

str_split(sigEnGS[1,]$EnrichmentDirection_ClusterID,"-")[[1]][2]

mastertable <- NULL
for (i in 1:length(sigEnGS[,1])){
  if(i == 119) next
  print(i)
  print(list.files( pattern = glob2rx(paste(sigEnGS[i,]$Comparison,"_","Cluster_",str_split(sigEnGS[i,]$EnrichmentDirection_ClusterID,"-")[[1]][2],"_ExpMatrix_Calvin_manual_genesets.",sigEnGS[i,]$GS,".report.",sigEnGS[i,]$EnrichmentDirection_ClusterID,"*.txt",sep = "")) ))
  filetxt <- list.files( pattern = glob2rx(paste(sigEnGS[i,]$Comparison,"_","Cluster_",str_split(sigEnGS[i,]$EnrichmentDirection_ClusterID,"-")[[1]][2],"_ExpMatrix_Calvin_manual_genesets.",sigEnGS[i,]$GS,".report.",sigEnGS[i,]$EnrichmentDirection_ClusterID,"*.txt",sep = "")) )
  table <- as.tibble(read.delim(filetxt,header = T)) %>% dplyr::filter(CORE_ENRICHMENT == "YES") %>% select(GENE,S2N) %>% add_column(Comparison=sigEnGS[i,]$Comparison, Cluster=str_split(sigEnGS[i,]$EnrichmentDirection_ClusterID,"-")[[1]][2],GeneSetName=sigEnGS[i,]$GS)
  mastertable <- rbind(mastertable,table)
  print(head(table))
  }

mn <- FindVariableGenes(mn, do.plot = F, display.progress = F)
hv.genes <- head(rownames(mn@hvg.info), 1000)
mn <- RunTSNE(mn, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)

TSNEPlot(mn, do.label = TRUE)
FeaturePlot(mn, features.plot = c("Pde4b"),cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
