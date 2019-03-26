####################################################################### 
#######    Imputation for Ultra-high-throughput scRNA-seq data
####### Requirement: conda 4.5.4, R 3.4.3, Rtsne 0.13, Rcpp 0.12.18, 
#######              Ramgic 1.4.0, scImpute 0.0.9
#######################################################################
source('function.R')

### Karaiskos
ImputationAndDimReduction(Raw_fileName='Karaiskos.Rdata',Cluster_fileName='Karaiskos_clust.Rdata',
	file='Karaiskos',Kcluster=4,cores=20,seed.value=12345)
plotFig(file='Karaiskos',Cluster_fileName='Karaiskos_clust.Rdata')

### PBMC
ImputationAndDimReduction(Raw_fileName='PBMC.Rdata',Cluster_fileName='PBMC_clust.Rdata',
	file='PBMC',Kcluster=1,cores=20,seed.value=12345)
plotFig(file='PBMC',Cluster_fileName='PBMC_clust.Rdata')

### Hrvatin
ImputationAndDimReduction(Raw_fileName='Hrvatin.Rdata',Cluster_fileName='Hrvatin_clust.Rdata',
	file='Hrvatin',,Kcluster=1,cores=20,seed.value=16431)
plotFig(file='Hrvatin',Cluster_fileName='Hrvatin_clust.Rdata')

### Alles
ImputationAndDimReduction(Raw_fileName='Alles.Rdata',Cluster_fileName='Alles_clust.Rdata',
	file='Alles',Kcluster=1,cores=20,seed.value=3454)
plotFig2(file='Alles',Cluster_fileName='Alles_clust.Rdata')
