## scRMD: Imputation for single cell RNA-seq data via restricted robust matrix decomposition
Chong Chen, Changjing Wu 2017-9-14

## Introduction
scRMD is developed to impute single cell RNA data with dropouts. scRMD assumes the the underlying
expression profile of genes is low rank and the dropout events are rare compared with true zero expression.

## Installation
scRMD can be installed by simplely run:

``` r
install.packages("devtools")         
library(devtools)           
install_github("ruibinxx/scRMD")
```

## Quick start

``` r
library(scRMD)
set.seed(2017)
K=3; Kn=50; Ndiff=100; Nsame=10000; logMean=1.8; logSd=0.5; 
ZeroRate = 0.5; sigmahetero = 0.1; sigmahomo = 0.2; drbase = 1; dr = 0.2;
sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "cluster")
cutoff = quantile(sData$de[sData$de>0], 0.05)
res.rmd <- rmd(sData$de, candidate = cutoff)
pca.rmd <- prcomp(res.rmd$exprs)
cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
CalculateARI(sData$label, cl.rmd$cluster)
```
