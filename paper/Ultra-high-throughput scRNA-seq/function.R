### functions
GeneNorm <- function(Y, percent = 0.05) {
  # Y is the observed p by n raw counts, p is #gene, n is #cell
  # Genes expressed in less than percent cells are excluded
  # The raw counts are normalized by CPM
  # Return the log10(counts+1) matrix
  n <- dim(Y)[2]
  gene.exprs.count <- rowSums(Y != 0)
  Y = Y[gene.exprs.count > n * percent, ]
  Y = t(t(Y)/rowSums(t(Y))*1000000)
  return (log10(Y+1))
}

DR <- function(X,Type){
  # X is n by k dimension resuction result, n is #cell
  # Tye is the true clusters for each cell
  IN <- OUT <- c()
  D = as.matrix(dist(X))
  n = dim(D)[1]
    U = unique(Type)
    for(u in U){
  index = which(Type == u)
  IN = c(IN,as.vector(D[index,index]))
  OUT = c(OUT,as.vector(D[index,-index]))
    }
  d = c(mean(IN[IN!=0])/mean(OUT),median(IN[IN!=0])/median(OUT))
  names(d) = c('Ratio of mean','Ratio of median')
  return(d)
}

### Imputation function
ImputationAndDimReduction <- function(Raw_fileName,Cluster_fileName,file,Location = getwd(),d=50,
  Kcluster=1,cores=20,seed.value=12345){
  # Location: Location of data
  # Raw_fileName: file name of raw data, a Rdata file which recorded UMI counts
  # Cluster_fileName: file name of cluster data
  # file: Output file name
  # Kcluster and cores are paramters for scImpute
  # seed.value: the seed value of Rtsne

  ## packages and functions 
  require(scRMD)
  require(Rtsne)
  require(Rmagic)
  require(scImpute)
  setwd(Location)

  ## load data
  load(Raw_fileName)
  load(Cluster_fileName)

  ## Normalize data
  mat = as.matrix(mat)
  res.raw = GeneNorm(mat)
  res.raw = t(res.raw)

  ## Imputation
  # scRMD
  cutoff = quantile(res.raw[res.raw!=0],0.05)
  res.rmd = rmd(res.raw,candidate = cutoff)
  save(res.rmd,file=paste0(file,'_rmd.Rdata'))
  print('Finish scRMD')
  # MAGIC
  res.mag = magic(res.raw)
  save(res.mag,file=paste0(file,'_mag.Rdata'))
  print('Finish MAGIC')
  # scImpute
  write.csv(t(10^res.raw+1),file='RawCount.csv')
  scimpute(# full path to raw count matrix
         count_path = 'RawCount.csv', 
         infile = 'csv',           # format of input file
         outfile = 'csv',          # format of output file
         out_dir = './',           # full path to output directory
         labeled = FALSE,          # cell type labels not available
         drop_thre = 0.5,          # threshold set on dropout probability
         Kcluster = Kcluster,             # 2 cell subpopulations
         ncores = cores) 
  res.sci = read.csv('scimpute_count.csv')
  res.sci = as.matrix(res.sci[,-1])
  res.sci = t(GeneNorm(res.sci))
  print('Finish scImpute')

  ## Percentage of zero factors
  n = dim(res.raw)[1]
  p = dim(res.raw)[2]
  raw.zero = sum(res.raw == 0)/n/p
  rmd.zero = sum(res.rmd$exprs == 0)/n/p
  rmd.S.zero = sum(res.rmd$S == 0)/n/p
  mag.zero = sum(res.mag$result == 0)/n/p
  sci.zero = sum(res.sci == 0)/n/p
  p0 = c(raw.zero,rmd.zero,mag.zero,sci.zero,rmd.S.zero)
  names(p0) = c('Raw','scRMD','MAGIC','scImpute','scRMD(S)')
  write.table(p0,file=paste0(file,'_p0.txt'),quote=FALSE)
  print('Finish caculating the percentage of zero factors')

  ###################################################################################
  ## tSNE
  set.seed(seed.value)
  raw.tsne = Rtsne(res.raw,initial_dims = d)$Y
  set.seed(seed.value)
  rmd.tsne = Rtsne(res.rmd$exprs,initial_dims = d)$Y
  set.seed(seed.value)
  mag.tsne = Rtsne(res.mag$result,initial_dims = d)$Y
  set.seed(seed.value)
  sci.tsne = Rtsne(res.sci,initial_dims = d)$Y
  result = list(raw.tsne,rmd.tsne,mag.tsne,sci.tsne)
  names(result) = c('Raw','scRMD','MAGIC','scImpute')
  save(result,file=paste0(file,'_tsne.Rdata'))
  # Ratio of between and within clusters
  Type = clust$Cluster
  index = which(Type != 'others')
  raw.ratio = DR(raw.tsne[index,],Type[index])
  rmd.ratio = DR(rmd.tsne[index,],Type[index])
  mag.ratio = DR(mag.tsne[index,],Type[index])
  sci.ratio = DR(sci.tsne[index,],Type[index])
  Ratio = rbind(raw.ratio,rmd.ratio,mag.ratio,sci.ratio)
  colnames(Ratio) =  c('mean','median')
  rownames(Ratio) = c('Raw','scRMD','MAGIC','scImpute')
  write.table(Ratio,file=paste0(file,'_tsne_Ratio.txt'),quote=FALSE)
  print('Finish tSNE')
}


plotFig <- function(file,Cluster_fileName,Location=getwd()){
  ## Color and legend
  setwd(Location)  
  clust = clust[index,]
  leg.col = c()
  leg.txt = sort(unique(as.character(clust$Cluster)))
  for(u in leg.txt){
    tmp = clust$Color[clust$Cluster == u]
    leg.col = c(leg.col,as.character(tmp[1]))
  }
  COLOR = as.character(clust$Color)
  
  ### Type of points
  if(length(clust$Pch) ==0){
    PCH = rep(20,length(COLOR))
    leg.pch = rep(20,length(leg.txt))
  }else{
      leg.pch = c()
      leg.txt = sort(unique(as.character(clust$Cluster)))
      for(u in leg.txt){
        tmp = clust$Pch[clust$Cluster == u]
        leg.pch = c(leg.pch,as.numeric(tmp))}
  }

  ### plot
  load(file=paste0(file,'_tsne.Rdata'))
  Type = c('Raw','scRMD','MAGIC','scImpute')
  for(i in 1 : 4){
  png(file = paste0(file,'_',Type[i],'.png'),width=800,height=800)
  plot(result[[i]],col=COLOR,pch=PCH,xlab='',ylab='',cex.axis=2,main=Type[i])
  dev.off()
  }
  png(file = paste0(file,'_Legend.png'),width=800,height=800)
  plot(c(0,20),c(0,20),type='n',xlab='',ylab='',bty='n',yaxt='n',xaxt='n')
  legend(c(0,20),legend=leg.txt,col=leg.col,pch = leg.pch,bty='n')
  dev.off()
}
