#Please change these according to your paths.
.libPaths(c(.libPaths(),"/n/home13/yasinkaymaz/biotools/Rlibs/"))
workingDirectory='/n/home13/yasinkaymaz/LabSpace/testdata/DulacLab/GSEA/MDT_mn'
source("~/codes/Brainformatics/scripts/functions.R")

### Run
library(Matrix)
library(Seurat)
library(dplyr)
library(tidyverse)
load(paste(workingDirectory,"mn.RObj",sep=""))
#load(paste(workingDirectory,"ACC.Neurons.RObj",sep=""))

setwd(workingDirectory)

setwd("~/Documents/DulacLab/ACC_an/")


load("an.RObj")
source("~/Dropbox/codes/Brainformatics/scripts/functions.R")
setwd("~/Documents/DulacLab/MDT_mn/")
load("mn.RObj")
mmDatasets=c("Calvin_manual_genesets.gmt")
i=1
exp.Seu.obj=mn
for (ds in mmDatasets){
  
  for (exp.Seu.obj in c(mn)){
    i=i+1
    #RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Subordinate", GeneSet = ds, outputDir = './')
    #RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Control", GeneSet = ds, outputDir = './')
    #RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Subordinate", Cond2 = "Control", GeneSet = ds, outputDir = './')
    
    Sig.Enrichment.ESTable <- SummarizeGSEAoutputs(GSEAoutputDir = "~/Documents/DulacLab/MDT_mn/")
      
    DEGs.dc <- RunDGEA(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Control")
    DEGs.sd <- RunDGEA(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Subordinate")
    DEGs.sc <- RunDGEA(SeuratObj = exp.Seu.obj, Cond1 = "Subordinate", Cond2 = "Control")
    
    #Combine results:
    DEGs <- dplyr::bind_rows(DEGs.dc, DEGs.sd, DEGs.sc, .id='id')
    DEGs$fdr <- p.adjust(DEGs$p_val, method = 'fdr')
    DEGs$id <- factor(DEGs$id)
    levels(DEGs$id) <- list('Dominant_Control'='1', 'Dominant_Subordinate'='2', 'Subordinate_Control'='3')
    DEGs <- arrange(DEGs,fdr)
    save(DEGs,file = paste("DEGs",i,"Rdata",sep="."))
    merged <- inner_join(Sig.Enrichment.ESTable, DEGs, by = c("GENE" = "gene", "Comparison" = "id", "cluster" = "cluster" )) %>% as.data.frame()
    write_tsv(merged,file.path(paste("GSEA-DGE-mergedtables.significantEnrichments",i,"txt",sep=".")))
  }
  
}

identical(table,merged)

table <- read.delim("GSEA-DGE-mergedtables.significantEnrichments.8.txt",header = T)
library(plotly)
library(gridExtra)
library(htmlwidgets)


comparisons= c("Dominant_Subordinate","Dominant_Control","Subordinate_Control")

for (comp in comparisons){
  p <- table %>% as.tibble() %>% filter(Comparison == comp) %>% 
    ggplot(aes(x=-log10(fdr), y=S2N, color=as.factor(cluster) )) + 
    geom_point(aes(text=sprintf("Gene: %s<br>GeneSet: %s<br>Cluster: %s", GENE, GS, cluster))) + 
    geom_vline(xintercept = -log10(0.05), linetype="dotted", color = "blue", size=0.5)+ 
    geom_hline(yintercept = 0, linetype="dotted", color = "red", size=0.5)+
    ggtitle(comp)
  pp <- ggplotly(p)
  htmlwidgets::saveWidget(widget=pp,paste(comp,".html",sep=""), selfcontained = FALSE)
}

comp="Dominant_Subordinate"
p <- table %>% as.tibble() %>% filter(Comparison == comp) %>% 
  ggplot(aes(x=-log10(p_val), y=S2N, color=as.factor(cluster) )) + 
  geom_point(aes(text=sprintf("Gene: %s<br>GeneSet: %s<br>Cluster: %s", GENE, GS, cluster))) + 
  geom_vline(xintercept = -log10(0.05), linetype="dotted", color = "blue", size=0.5)+ 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=0.5)+
  ggtitle(comp)

ggplotly(p)








hv.genes <- head(rownames(an@hvg.info), 1000)
an <- ScaleData(an, 
                genes.use = hv.genes, 
                vars.to.regress = c("nUMI", "percent.mito"),
                display.progress = FALSE)
an <- RunPCA(an, 
             pc.genes = hv.genes, 
             do.print = FALSE)

an <- JackStraw(an, num.pc = 60, display.progress = T, do.par = T, num.cores = 6, maxit = 1000)

JackStrawPlot(an, PCs=1:20)
PCElbowPlot(an, num.pc = 20)

an <- RunTSNE(an, dims.use = 1:20, do.fast = TRUE)

pdf("tsne.pdf",width = 10,height = 10)
TSNEPlot(an, do.label = TRUE, group.by="tree.ident")
FeaturePlot(an, features.plot = c("Snhg11"))
dev.off()







load("~/Documents/DulacLab/Eric-Adam/an.RObj")

source("~/Dropbox/codes/Brainformatics/scripts/functions.R")

DEGs.sd <- RunDGEA(SeuratObj = an, Cond1 = "Dominant", Cond2 = "Subordinate") 
library(tidyverse)
DEGs.sd %>% as.tibble() %>% filter(p_val_adj < 0.05) %>% arrange(cluster)

DEGs.sd.b <- RunDGEA(SeuratObj = an, Cond1 = "Dominant", Cond2 = "Subordinate") 
DEGs.sd.b %>% as.tibble() %>% filter(p_val_adj < 0.05) %>% arrange(cluster)








GSEAoutputDir = "~/Documents/DulacLab/MDT_mn/"
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
      print(geneSetReportfile)
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


