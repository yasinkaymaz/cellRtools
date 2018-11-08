
library(tidyverse)

#Please change this directory   !!!
setwd("~/Downloads/GSEA_MDT/")

mySumFiles <- list.files(pattern="*SUMMARY.RESULTS.REPORT*")

majorSummaryTable <- NULL

for (i in 1:length(mySumFiles)){
  print(i)
  sumTable <- read.delim(mySumFiles[i]) %>% as.tibble() %>% add_column(Comparison=strsplit(mySumFiles[i],"_Clust")[[1]][1],EnrichmentDirection_ClusterID=strsplit(mySumFiles[i],"\\.")[[1]][5])
  majorSummaryTable <- bind_rows(majorSummaryTable, sumTable)
}

majorSummaryTable <- majorSummaryTable %>% as.tibble() %>% mutate(pAdj.Nom=p.adjust(NOM.p.val,method="BH")) %>% arrange(pAdj.Nom)

majorSummaryTable %>% write_tsv(file.path("All.Enrichment.stats.txt"))

sigtable <- majorSummaryTable %>% filter(pAdj.Nom < 0.05) %>% unite(plotfilename, Comparison,GS,EnrichmentDirection_ClusterID,sep="*",remove = FALSE)

majorSummaryTable %>% filter(pAdj.Nom < 0.05) %>% unite(plotfilename, Comparison,GS,EnrichmentDirection_ClusterID,sep="*",remove = FALSE) %>% write_tsv(file.path("significant.Enrichments.txt"))

#Analysis plots

pdf("Summary.Plots.pdf")
sigtable <- sigtable %>% separate(EnrichmentDirection_ClusterID,into=c("EnrichmentDirection","ClusterID"),sep="-")

ggplot(sigtable, aes(ClusterID,fill=EnrichmentDirection))+geom_histogram(stat="count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(sigtable, aes(Comparison, fill=EnrichmentDirection))+geom_histogram(stat="count",position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

mdt.clusters = c(6,15,17,19,20,22,25) 

sigtable.mdt.clusters <- sigtable %>% filter(ClusterID %in% mdt.clusters)

ggplot(sigtable.mdt.clusters, aes(ClusterID,fill=EnrichmentDirection))+geom_histogram(stat="count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(sigtable.mdt.clusters, aes(Comparison, fill=EnrichmentDirection))+geom_histogram(stat="count",position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

pdf("GS.counts.pdf",width = 10,height = 10)

ggplot(sigtable,aes(GS,fill=ClusterID))+geom_histogram(stat="count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(sigtable.mdt.clusters,aes(GS,fill=ClusterID))+geom_histogram(stat="count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
