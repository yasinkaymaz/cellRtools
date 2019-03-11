### Cell Type Prediction

The goal of this project is to be able to determine major cell types of cell samples in the single cell RNAseq (scRNAseq) datasets. Common methods for deciding type of the cells often involve manually checking known tissue/cell specific gene expression levels. Sensitive methods for automatically determining cell types is currently lacking. Therefore, our effort is to develop a machine learning approach to predict cell labels using gene expression levels from scRNAseq datasets. This will allow researchers to find out existing cell types in their experimental outcomes and use these predictions to further fine tune their data for downstream analysis, such as marker gene selection, differential expression test, etc.


#### Training Data Preparation
1. PCA
2. Gene selection


####

# Hierarchical Random Forest

1.  Build a Hierarchical clustering / Phylogenetic tree
2.  Create Local Classifiers at each node
    1.  Split the data into node-subclusters
3.  Build LC Map. Store them into an easily accessible format.
4.  Build Probability functions for each leaf


#### [Example Tutorial](tutorial_1.md)
