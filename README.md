# cellRtools
This repository is a collection of analytical methods for single cell sequencing data.

### Install

To be able to use the functions, simply clone the repository:
```{bash}
git clone https://github.com/yasinkaymaz/cellRtools.git
```

### Overview

#### [Cell Type Prediction](Docs/CellTyper.md)

The goal of this project is to be able to determine major cell types of cell samples in the single cell RNAseq (scRNAseq) datasets. Common methods for deciding type of the cells often involve manually checking known tissue/cell specific gene expression levels. Sensitive methods for automatically determining cell types is currently lacking. Therefore, our effort is to develop a machine learning approach to predict cell labels using gene expression levels from scRNAseq datasets. This will allow researchers to find out existing cell types in their experimental outcomes and use these predictions to further fine tune their data for downstream analysis, such as marker gene selection, differential expression test, etc.


#### Gene Set Enrichment Analysis

This project is for providing a user friendly method to perform Gene Set Enrichment Analysis between clusters (or subsets of cells in clusters). The main GSEA software has been adapted from [Broad Institute](http://software.broadinstitute.org/gsea/) and restored for compatibility with latest R distributions.
