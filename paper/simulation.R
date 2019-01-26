library("kernlab")
library("Rmagic")
library("reticulate")
library("scImpute")
library("data.table")
source("Scimpute/simulation.R")
source("Scimpute/calculate_weight.R")
source("Scimpute/dmix.R")
source("Scimpute/get_mix_parameters.R")
source("Scimpute/imputation_model.R")
source("Scimpute/rmix.R")
source("function.R")
## Tables


use_python("/data1/XiLab/chenchong/anaconda3/bin/python")
################################################ Clustering Analysis ##################################################

K=3; Kn=50; Ndiff = 100; Nsame=10000; logMean=1.8; logSd=0.5; nsium = 50; num.cores = 10;
candidate_thres <- c(0, 0.02, 0.05, 0.08, 0.1)
## Tables
## Changing sigmahomo
set.seed(2017)
ZeroRate = 0.3; sigmahetero = 0.1; sigmahomo = c(0,0.1,0.2,0.3,0.4); drbase = 1; dr = 0.2;
simu.cluster <- data.frame() 
for (i in 1:5){
	for(k in 1:nsium){
		print(k)
    # df_tmp = data.frame(method = rep(c("Oracle","RAW","RMD","scImpute","MAGIC")), 
    #   SigmaHomo= rep(sigmahomo[i], 5), ari = numeric(5))
		sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo[i], sigmahetero, type = "cluster")
    fwrite(data.table(t(sData$te)), paste("simulation_data/oracle_sigmahomo", sigmahomo[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    fwrite(data.table(t(sData$de)), paste("simulation_data/raw_sigmahomo", sigmahomo[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
  	res.sci = imputation_model8(count = t(sData$de+log10(1.01)), labeled = False, 
  			point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("simulation_data/sci_sigmahomo", sigmahomo[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
  	res.magic <- magic(sData$de, genes="all_genes")$result
    fwrite(data.table(t(res.magic)), paste("simulation_data/magic_sigmahomo", sigmahomo[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], qq))
      fwrite(data.table(t(res.rmd$exprs)), paste("simulation_data/rmd_sigmahomo", sigmahomo[i], "_rep", k, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    }
  # 	res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  # 	pca.oracle <- prcomp(sData$te)
  # 	pca.raw <- prcomp(sData$de)
		# pca.rmd <- prcomp(res.rmd$exprs)
		# pca.sci <- prcomp(res.sci)
		# pca.magic <- prcomp(res.magic)
		# cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
		# cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
		# cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
		# cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
		# cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
		# df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
		# df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
		# df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
		# df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
		# df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
		# simu.cluster = rbind(simu.cluster, df_tmp)
	}
}
#write.csv(simu.cluster,file="Result/sigmahomo.csv",quote=F,row.names=F,col.names=T)




## Tables
## Changing sigmahetero

set.seed(2017)
ZeroRate = 0.5; sigmahetero = c(0,0.05,0.1,0.15,0.2); sigmahomo = 0.2; drbase = 1; dr = 0.2;
simu.cluster <- data.frame()
for (i in 1:5){
	for (k in 1:nsium){
		# df_tmp = data.frame(method = rep(c("Oracle","RAW","RMD","scImpute","MAGIC")), 
		# 	SigmaHetero= rep(sigmahetero[i], 5), ari = numeric(5))
		sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero[i], type = "cluster")
    fwrite(data.table(t(sData$te)), paste("simulation_data/oracle_sigmahetro", sigmahetro[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    fwrite(data.table(t(sData$de)), paste("simulation_data/raw_sigmahetro", sigmahetro[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    res.sci = imputation_model8(count = t(sData$de+log10(1.01)), labeled = False, 
        point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("simulation_data/sci_sigmahetro", sigmahetro[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    res.magic <- magic(sData$de, genes="all_genes")$result
    fwrite(data.table(t(res.magic)), paste("simulation_data/magic_sigmahetro", sigmahetro[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], qq))
      fwrite(data.table(t(res.rmd$exprs)), paste("simulation_data/rmd_sigmahetro", sigmahetro[i], "_rep", k, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    }
  # 	res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  # 	pca.oracle <- prcomp(sData$te)
  # 	pca.raw <- prcomp(sData$de)
		# pca.rmd <- prcomp(res.rmd$exprs)
		# pca.sci <- prcomp(res.sci)
		# pca.magic <- prcomp(res.magic)
		# cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
		# cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
		# cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
		# cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
		# cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
		# df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
		# df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
		# df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
		# df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
		# df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
		# simu.cluster = rbind(simu.cluster, df_tmp)

	}
}
#write.csv(simu.cluster,file="Result/sigmahetero.csv",quote=F,row.names=F,col.names=T)


## Tables
## Changing ZeroRate

set.seed(2017)
ZeroRate = c(0.3,0.4,0.5,0.6,0.7); sigmahetero = 0.1; sigmahomo = 0.2; drbase = 1; dr = 0.2;
simu.cluster <- data.frame()
for (i in 1:8){
	for (k in 1:nsium){
    # df_tmp = data.frame(method = rep(c("Oracle","RAW","RMD","scImpute","MAGIC")), 
    #   zerorate= rep(ZeroRate[i], 5), ari = numeric(5))
		sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate[i], drbase, dr, sigmahomo, sigmahetero, type = "cluster")
    fwrite(data.table(t(sData$te)), paste("simulation_data/oracle_ZeroRate", ZeroRate[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    fwrite(data.table(t(sData$de)), paste("simulation_data/raw_ZeroRate", ZeroRate[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    res.sci = imputation_model8(count = t(sData$de+log10(1.01)), labeled = False, 
        point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("simulation_data/sci_ZeroRate", ZeroRate[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    res.magic <- magic(sData$de, genes="all_genes")$result
    fwrite(data.table(t(res.magic)), paste("simulation_data/magic_ZeroRate", ZeroRate[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], qq))
      fwrite(data.table(t(res.rmd$exprs)), paste("simulation_data/rmd_ZeroRate", ZeroRate[i], "_rep", k, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    }
 #  	res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
 #  	pca.oracle <- prcomp(sData$te)
 #  	pca.raw <- prcomp(sData$de)
	# 	pca.rmd <- prcomp(res.rmd$exprs)
	# 	pca.sci <- prcomp(res.sci)
	# 	pca.magic <- prcomp(res.magic)
	# 	cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
	# 	cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
	# 	cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
	# 	cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
	# 	cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
	# 	df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
	# 	df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
	# 	df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
	# 	df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
	# 	df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
	# 	simu.cluster = rbind(simu.cluster, df_tmp)
	}
}
#write.csv(simu.cluster,file="Result/ZeroRate.csv",quote=F,row.names=F,col.names=T)



## Tables
## Changing dr

set.seed(2017)
ZeroRate = 0.5; sigmahetero = 0.1; sigmahomo = 0.2; drbase = 1; dr = c(0.1,0.2,0.3,0.4,0.5);
simu.cluster <- data.frame()
for (i in 1:5){
	for (k in 1:nsium){
    # df_tmp = data.frame(method = rep(c("Oracle","RAW","RMD","scImpute","MAGIC")), 
    #   DR= rep(dr[i], 5), ari = numeric(5))
		sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr[i], sigmahomo, sigmahetero, type = "cluster")
    fwrite(data.table(t(sData$te)), paste("simulation_data/oracle_dr", dr[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    fwrite(data.table(t(sData$de)), paste("simulation_data/raw_dr", dr[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    res.sci = imputation_model8(count = t(sData$de+log10(1.01)), labeled = False, 
        point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("simulation_data/sci_dr", dr[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    res.magic <- magic(sData$de, genes="all_genes")$result
    fwrite(data.table(t(res.magic)), paste("simulation_data/magic_dr", dr[i], "_rep", k, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], qq))
      fwrite(data.table(t(res.rmd$exprs)), paste("simulation_data/rmd_dr", dr[i], "_rep", k, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    }
  	# res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  # 	pca.oracle <- prcomp(sData$te)
  # 	pca.raw <- prcomp(sData$de)
		# pca.rmd <- prcomp(res.rmd$exprs)
		# pca.sci <- prcomp(res.sci)
		# pca.magic <- prcomp(res.magic)
		# cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
		# cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
		# cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
		# cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
		# cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
		# df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
		# df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
		# df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
		# df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
		# df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
		# simu.cluster = rbind(simu.cluster, df_tmp)
	}
}
# write.csv(simu.cluster,file="Result/dr.csv",quote=F,row.names=F,col.names=T)


######################################## sensitivity analysis ###############################################################################################

K=3; Kn=50; Ndiff = 100; Nsame=10000; logMean=1.8; logSd=0.5; nsium = 50; num.cores = 10;
candidate_thres <- c(0, 0.02, 0.05, 0.1, 0.2)
## Tables
## Changing sigmahomo
set.seed(2017)
ZeroRate = 0.3; sigmahetero = 0.1; sigmahomo = c(0,0.1,0.2,0.3,0.4); drbase = 1; dr = 0.2;
simu.cluster <- data.frame() 
for (i in 1:5){
  for(k in 1:nsium){
    print(k)
    df_tmp = data.frame(method = rep(c("cutoff_0","cutoff_0.02","cutoff_0.05","cutoff_0.1","cutoff_0.2")), 
      SigmaHomo= rep(sigmahomo[i], 5), ari = numeric(5))
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo[i], sigmahetero, type = "cluster")
    for (j in c(1:5)){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], candidate_thres[j]))
      pca.rmd <- prcomp(res.rmd$exprs)
      cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
      df_tmp$ari[j]  = ARI(sData$label, cl.rmd$cluster)
    }
  #   res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  #   pca.oracle <- prcomp(sData$te)
  #   pca.raw <- prcomp(sData$de)
    # pca.rmd <- prcomp(res.rmd$exprs)
    # pca.sci <- prcomp(res.sci)
    # pca.magic <- prcomp(res.magic)
    # cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
    # cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
    # cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
    # cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
    # cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
    # df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
    # df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
    # df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
    # df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
    # df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
    simu.cluster = rbind(simu.cluster, df_tmp)
  }
}
write.csv(simu.cluster,file="Result/sigmahomo_sensitivity.csv",quote=F,row.names=F,col.names=T)




## Tables
## Changing sigmahetero

set.seed(2017)
ZeroRate = 0.5; sigmahetero = c(0,0.05,0.1,0.15,0.2); sigmahomo = 0.2; drbase = 1; dr = 0.2;
simu.cluster <- data.frame()
for (i in 1:5){
  for (k in 1:nsium){
    df_tmp = data.frame(method = rep(c("cutoff_0","cutoff_0.02","cutoff_0.05","cutoff_0.1","cutoff_0.2")), 
      SigmaHetero= rep(sigmahetero[i], 5), ari = numeric(5))
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero[i], type = "cluster")
      for (j in c(1:5)){
          res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], candidate_thres[j]))
        pca.rmd <- prcomp(res.rmd$exprs)
          cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
          df_tmp$ari[j]  = ARI(sData$label, cl.rmd$cluster)
      }
  #   res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  #   pca.oracle <- prcomp(sData$te)
  #   pca.raw <- prcomp(sData$de)
    # pca.rmd <- prcomp(res.rmd$exprs)
    # pca.sci <- prcomp(res.sci)
    # pca.magic <- prcomp(res.magic)
    # cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
    # cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
    # cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
    # cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
    # cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
    # df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
    # df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
    # df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
    # df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
    # df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
    simu.cluster = rbind(simu.cluster, df_tmp)
  }
}
write.csv(simu.cluster,file="Result/sigmahetero_sensitivity.csv",quote=F,row.names=F,col.names=T)


## Tables
## Changing ZeroRate

set.seed(2017)
ZeroRate = c(0.3,0.4,0.5,0.6,0.7); sigmahetero = 0.1; sigmahomo = 0.2; drbase = 1; dr = 0.2;
simu.cluster <- data.frame()
for (i in 1:8){
  for (k in 1:nsium){
    df_tmp = data.frame(method = rep(c("cutoff_0","cutoff_0.02","cutoff_0.05","cutoff_0.1","cutoff_0.2")), 
      zerorate= rep(ZeroRate[i], 5), ari = numeric(5))
  sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate[i], drbase, dr, sigmahomo, sigmahetero, type = "cluster")
    for (j in c(1:5)){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], candidate_thres[j]))
      pca.rmd <- prcomp(res.rmd$exprs)
      cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
      df_tmp$ari[j]  = ARI(sData$label, cl.rmd$cluster)
    }
  #   res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  #   pca.oracle <- prcomp(sData$te)
  #   pca.raw <- prcomp(sData$de)
    # pca.rmd <- prcomp(res.rmd$exprs)
    # pca.sci <- prcomp(res.sci)
    # pca.magic <- prcomp(res.magic)
    # cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
    # cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
    # cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
    # cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
    # cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
    # df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
    # df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
    # df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
    # df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
    # df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
    simu.cluster = rbind(simu.cluster, df_tmp)
  }
}
write.csv(simu.cluster,file="Result/ZeroRate_sensitivity.csv",quote=F,row.names=F,col.names=T)



## Tables
## Changing dr

set.seed(2017)
ZeroRate = 0.5; sigmahetero = 0.1; sigmahomo = 0.2; drbase = 1; dr = c(0.1,0.2,0.3,0.4,0.5);
simu.cluster <- data.frame()
for (i in 1:5){
  for (k in 1:nsium){
    df_tmp = data.frame(method = rep(c("cutoff_0","cutoff_0.02","cutoff_0.05","cutoff_0.1","cutoff_0.2")), 
      DR= rep(dr[i], 5), ari = numeric(5))
  sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr[i], sigmahomo, sigmahetero, type = "cluster")
    for (j in c(1:5)){
      res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], candidate_thres[j]))
      pca.rmd <- prcomp(res.rmd$exprs)
      cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
      df_tmp$ari[j]  = ARI(sData$label, cl.rmd$cluster)
    }
  #   res.rmd <- rmd(sData$de, candidate = quantile(sData$de[sData$de > 0], 0.05))
  #   pca.oracle <- prcomp(sData$te)
  #   pca.raw <- prcomp(sData$de)
    # pca.rmd <- prcomp(res.rmd$exprs)
    # pca.sci <- prcomp(res.sci)
    # pca.magic <- prcomp(res.magic)
    # cl.oracle <- kmeans(pca.oracle$x[,1:2],K,nstart = 100)
    # cl.raw <- kmeans(pca.raw$x[,1:2],K,nstart = 100)
    # cl.rmd <- kmeans(pca.rmd$x[,1:2],K,nstart = 100)
    # cl.sci <- kmeans(pca.sci$x[,1:2],K,nstart = 100)
    # cl.magic <- kmeans(pca.magic$x[,1:2],K,nstart = 100)
    # df_tmp$ari[1] = ARI(sData$label, cl.oracle$cluster)
    # df_tmp$ari[2]  = ARI(sData$label, cl.raw$cluster)
    # df_tmp$ari[3]  = ARI(sData$label, cl.rmd$cluster)
    # df_tmp$ari[4]  = ARI(sData$label, cl.sci$cluster)
    # df_tmp$ari[5]  = ARI(sData$label, cl.magic$cluster)
    simu.cluster = rbind(simu.cluster, df_tmp)
  }
}
write.csv(simu.cluster,file="Result/dr_sensitivity.csv",quote=F,row.names=F,col.names=T)



################################################### Differential Analysis ################################################
# Tables
# change dr, sigma, and Niff. dr=0.2,~70% 0; dr = 0.5,~50% 0. sigma = 0.3 or 0.1. Ndiff = 30 or 100
# record: method; zero prop; Fnorm; 0.01 pre,recal,Fscore,MCC; 0.05 same; time
# dr = 0.2, sigma = 0.1, Ndiff = 30
{
set.seed(2017)
K=2; Kn=50; Ndiff = 30; Nsame = 10000; logMean = 2; logSd = 1; 


nsimu <- 50
simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                      Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                      Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                      recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                      time = numeric(nsimu*5), stringsAsFactors = F)
for (k in 1:nsimu) {
  sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
  
  tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
  tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
  simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
  
  simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
  simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
  simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
  simu.df$method[(k-1)*5+4] <- "SCI"
  t1 <- proc.time()
  plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
  res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
  res.sci = t(res.sci - log10(1.01))
  simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
  simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
  
  simu.df$Fnorm[(k-1)*5+1] <- 0
  simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
  simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
  simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
  simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
  
  te.t.test <- c()
  for (i in 1:sData$P) {
    te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
  }
  te.p.adj <- p.adjust(te.t.test, method = "BH")
  de.t.test <- c()
  for (i in 1:sData$P) {
    de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
  }
  de.p.adj <- p.adjust(de.t.test, method = "BH")
  rmd.t.test <- c()
  for (i in 1:sData$P) {
    rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
  }
  rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
  sci.t.test <- c()
  for (i in 1:sData$P) {
    sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
  }
  sci.p.adj <- p.adjust(sci.t.test, method = "BH")
  magic.t.test <- c()
  for (i in 1:sData$P) {
    magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
  }
  magic.p.adj <- p.adjust(magic.t.test, method = "BH")
  
  golden <- c(rep(0,Ndiff), rep(1, Nsame))
  te.test.eval1 <- test.eval(golden, te.p.adj)
  de.test.eval1 <- test.eval(golden, de.p.adj)
  rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
  sci.test.eval1 <- test.eval(golden, sci.p.adj)
  magic.test.eval1 <- test.eval(golden, magic.p.adj)
  te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
  de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
  rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
  sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
  magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
  
  simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
  simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
  simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
  simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
  simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
  
  simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
  simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
  simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
  simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
  simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
  
  simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
  simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
  simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
  simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
  simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
  
  simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
  simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
  simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
  simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
  simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
  
  simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
  simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
  simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
  simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
  simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore

  simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
  simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
  simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
  simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
  simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
  
  simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
  simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
  simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
  simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
  simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
  
  simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
  simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
  simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
  simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
  simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
  
  print(k)
}
write.csv(simu.df, "simu-dr0.2sigma0.1Ndiff30.csv")
}

# dr = 0.2, sigma = 0.3, Ndiff = 30
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 30; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.3; sigmahomo = 0.3; drbase = 1; dr = 0.2;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.2sigma0.3Ndiff30.csv")
}

# dr = 0.2, sigma = 0.3, Ndiff = 100
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 100; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.3; sigmahomo = 0.3; drbase = 1; dr = 0.2;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.2sigma0.3Ndiff100.csv")
}

# dr = 0.2, sigma = 0.1, Ndiff = 100
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 100; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.1; sigmahomo = 0.1; drbase = 1; dr = 0.2;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.2sigma0.1Ndiff100.csv")
}

# dr = 0.5, sigma = 0.1, Ndiff = 30
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 30; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.1; sigmahomo = 0.1; drbase = 1; dr = 0.5;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.5sigma0.1Ndiff30.csv")
}

# dr = 0.5, sigma = 0.3, Ndiff = 30
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 30; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.3; sigmahomo = 0.3; drbase = 1; dr = 0.5;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.5sigma0.3Ndiff30.csv")
}

# dr = 0.5, sigma = 0.3, Ndiff = 100
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 100; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.3; sigmahomo = 0.3; drbase = 1; dr = 0.5;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.5sigma0.3Ndiff100.csv")
}

# dr = 0.5, sigma = 0.1, Ndiff = 100
{
  set.seed(2017)
  K=2; Kn=50; Ndiff = 100; Nsame = 10000; logMean = 2; logSd = 1; 
  ZeroRate = 0.1; sigmahetero = 0.1; sigmahomo = 0.1; drbase = 1; dr = 0.5;
  nsimu <- 50
  simu.df <- data.frame(method = character(nsimu*5), t0.prop = numeric(nsimu*5), d0.prop = numeric(nsimu*5),
                        Fnorm = numeric(nsimu*5), precision1 = numeric(nsimu*5), recall1 = numeric(nsimu*5), 
                        Fscore1 = numeric(nsimu*5), MCC1 = numeric(nsimu*5), precision5 = numeric(nsimu*5), 
                        recall5 = numeric(nsimu*5), Fscore5 = numeric(nsimu*5), MCC5 = numeric(nsimu*5), 
                        time = numeric(nsimu*5), stringsAsFactors = F)
  for (k in 1:nsimu) {
    sData = sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "DE")
    
    tmp1 <- sum(sData$de == 0) / prod(dim(sData$te))
    tmp2 <- sum(sData$te == 0) / prod(dim(sData$te))
    simu.df$t0.prop[((k-1)*5+1):(k*5)] <- tmp2; simu.df$d0.prop[((k-1)*5+1):(k*5)] <- tmp1;
    
    simu.df$method[(k-1)*5+1] <- "TE"; simu.df$time[(k-1)*5+1] <- 0
    simu.df$method[(k-1)*5+2] <- "DE"; simu.df$time[(k-1)*5+2] <- 0
    simu.df$method[(k-1)*5+3] <- "RMD"; simu.df$time[(k-1)*5+3] <- system.time(res.rmd <- rmd(sData$de))[3]
    simu.df$method[(k-1)*5+4] <- "SCI"
    t1 <- proc.time()
    plist = get_mix_parameters(t(sData$de) + log10(1.01)); 
    res.sci = imputation_model1(count = t(sData$de+log10(1.01)), point = log10(1.01), plist, drop_thre = 0.5);
    res.sci = t(res.sci - log10(1.01))
    simu.df$time[(k-1)*5+4] <- proc.time()[3] - t1[3]
    simu.df$method[(k-1)*5+5] <- "MAGIC"; simu.df$time[(k-1)*5+5] <- system.time(res.magic <- magic(sData$de))[3]
    
    simu.df$Fnorm[(k-1)*5+1] <- 0
    simu.df$Fnorm[(k-1)*5+2] <- norm(sData$te - sData$de, "F")
    simu.df$Fnorm[(k-1)*5+3] <- norm(sData$te - res.rmd$exprs, "F")
    simu.df$Fnorm[(k-1)*5+4] <- norm(sData$te - res.sci, "F")
    simu.df$Fnorm[(k-1)*5+5] <- norm(sData$te - res.magic, "F")
    
    te.t.test <- c()
    for (i in 1:sData$P) {
      te.t.test[i] <- my.t.test(sData$te[sData$label == 1, i], sData$te[sData$label == 2, i])
    }
    te.p.adj <- p.adjust(te.t.test, method = "BH")
    de.t.test <- c()
    for (i in 1:sData$P) {
      de.t.test[i] <- my.t.test(sData$de[sData$label == 1, i], sData$de[sData$label == 2, i])
    }
    de.p.adj <- p.adjust(de.t.test, method = "BH")
    rmd.t.test <- c()
    for (i in 1:sData$P) {
      rmd.t.test[i] <- my.t.test(res.rmd$exprs[sData$label == 1, i], res.rmd$exprs[sData$label == 2, i])
    }
    rmd.p.adj <- p.adjust(rmd.t.test, method = "BH")
    sci.t.test <- c()
    for (i in 1:sData$P) {
      sci.t.test[i] <- my.t.test(res.sci[sData$label == 1, i], res.sci[sData$label == 2, i])
    }
    sci.p.adj <- p.adjust(sci.t.test, method = "BH")
    magic.t.test <- c()
    for (i in 1:sData$P) {
      magic.t.test[i] <- my.t.test(res.magic[sData$label == 1, i], res.magic[sData$label == 2, i])
    }
    magic.p.adj <- p.adjust(magic.t.test, method = "BH")
    
    golden <- c(rep(0,Ndiff), rep(1, Nsame))
    te.test.eval1 <- test.eval(golden, te.p.adj)
    de.test.eval1 <- test.eval(golden, de.p.adj)
    rmd.test.eval1 <- test.eval(golden, rmd.p.adj)
    sci.test.eval1 <- test.eval(golden, sci.p.adj)
    magic.test.eval1 <- test.eval(golden, magic.p.adj)
    te.test.eval5 <- test.eval(golden, te.p.adj, 0.05, 0.05)
    de.test.eval5 <- test.eval(golden, de.p.adj, 0.05, 0.05)
    rmd.test.eval5 <- test.eval(golden, rmd.p.adj, 0.05, 0.05)
    sci.test.eval5 <- test.eval(golden, sci.p.adj, 0.05, 0.05)
    magic.test.eval5 <- test.eval(golden, magic.p.adj, 0.05, 0.05)
    
    simu.df$precision1[(k-1)*5+1] <- te.test.eval1$precision
    simu.df$precision1[(k-1)*5+2] <- de.test.eval1$precision
    simu.df$precision1[(k-1)*5+3] <- rmd.test.eval1$precision
    simu.df$precision1[(k-1)*5+4] <- sci.test.eval1$precision
    simu.df$precision1[(k-1)*5+5] <- magic.test.eval1$precision
    
    simu.df$precision5[(k-1)*5+1] <- te.test.eval5$precision
    simu.df$precision5[(k-1)*5+2] <- de.test.eval5$precision
    simu.df$precision5[(k-1)*5+3] <- rmd.test.eval5$precision
    simu.df$precision5[(k-1)*5+4] <- sci.test.eval5$precision
    simu.df$precision5[(k-1)*5+5] <- magic.test.eval5$precision
    
    simu.df$recall1[(k-1)*5+1] <- te.test.eval1$TPR
    simu.df$recall1[(k-1)*5+2] <- de.test.eval1$TPR
    simu.df$recall1[(k-1)*5+3] <- rmd.test.eval1$TPR
    simu.df$recall1[(k-1)*5+4] <- sci.test.eval1$TPR
    simu.df$recall1[(k-1)*5+5] <- magic.test.eval1$TPR
    
    simu.df$recall5[(k-1)*5+1] <- te.test.eval5$TPR
    simu.df$recall5[(k-1)*5+2] <- de.test.eval5$TPR
    simu.df$recall5[(k-1)*5+3] <- rmd.test.eval5$TPR
    simu.df$recall5[(k-1)*5+4] <- sci.test.eval5$TPR
    simu.df$recall5[(k-1)*5+5] <- magic.test.eval5$TPR
    
    simu.df$Fscore1[(k-1)*5+1] <- te.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+2] <- de.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+3] <- rmd.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+4] <- sci.test.eval1$Fscore
    simu.df$Fscore1[(k-1)*5+5] <- magic.test.eval1$Fscore
    
    simu.df$Fscore5[(k-1)*5+1] <- te.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+2] <- de.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+3] <- rmd.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+4] <- sci.test.eval5$Fscore
    simu.df$Fscore5[(k-1)*5+5] <- magic.test.eval5$Fscore
    
    simu.df$MCC1[(k-1)*5+1] <- te.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+2] <- de.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+3] <- rmd.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+4] <- sci.test.eval1$MCC
    simu.df$MCC1[(k-1)*5+5] <- magic.test.eval1$MCC
    
    simu.df$MCC5[(k-1)*5+1] <- te.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+2] <- de.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+3] <- rmd.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+4] <- sci.test.eval5$MCC
    simu.df$MCC5[(k-1)*5+5] <- magic.test.eval5$MCC
    
    print(k)
  }
  write.csv(simu.df, "simu-dr0.5sigma0.1Ndiff100.csv")
}


############################# Result Summary ######################################################




