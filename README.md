# scRMD: Imputation for single cell RNA-seq data via restricted robust matrix decomposition
Chong Chen, Changjing Wu 2017-9-14

# Introduction
scRMD is developed to impute single cell RNA data with dropouts. scRMD assumes the the underlying
expression profile of genes is low rank and the dropout events are rare compared with true zero expression.

# Installation
scRMD can be installed by simplely run:
>>install.packages("devtools")
           
>>library(devtools)
           
>>install_github("ChongC1990/scRMD")
