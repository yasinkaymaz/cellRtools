### Tutorial for Cell Type Prediction


Get the data [PBMC single cell RNAseq data with 2,700 cell](https://www.dropbox.com/s/kwd3kcxkmpzqg6w/pbmc3k_final.rds?dl=0). The final data is in Seurat object format and processed by following the [tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)
[Test set data is from 10X Genomics](http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz)

First of all load the functions
```{r}
source("cellRtools/main/functions.R")
#load the pbmc dataset
pbmc3k <- get(load("pbmc3k_final.rds"))

```

Training model with the curated dataset.
```{r}

pmbc.rf <- CellTyperTrainer2(ExpressionData = pbmc3k@data,
  CellLabels = pbmc3k@meta.data$ClusterNames_0.6,
  run.name = "pbmc.rf",
  do.splitTest = F,
  PCs = 20,
  improve = T)

```

Load a Test Set. Positive control
```{r}
exp.data <- Read10X(data.dir="~/Downloads/pbmc10K/")
pbmc10 <- SeuratWrapper(ExpData = exp.data, ProjectLabel = "pbmc10K",Normalize = T,scale.only.var = T,PCs = 20)

```
Predict the cell types in the test set
```{r}

pbmc10K <- CellTyper2(SeuratObject = pbmc10K, model = pmbc.rf, outputFilename = "pred.pbmc10K")

```