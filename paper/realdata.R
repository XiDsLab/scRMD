setwd("/data1/XiLab/chenchong/scRMD")
library("kernlab")
library("Rmagic")
library("reticulate")
library("scImpute")
library("data.table")
library("doParallel")
library("cidr")
library("SIMLR")
source("Scimpute/simulation.R")
source("Scimpute/calculate_weight.R")
source("Scimpute/dmix.R")
source("Scimpute/get_mix_parameters.R")
source("Scimpute/imputation_model.R")
source("Scimpute/rmix.R")
source("function.R")



################## Basic Functions ####
subsample <- function(A,j,prob){# A is p by n
  print(j)
  p <- ncol(A)
  sample.mtx <- rep(0,p)
  nonzero.idx <- A[j,] != 0
  sample.mtx[nonzero.idx] <- rbinom(sum(nonzero.idx), as.matrix(round(A[j,nonzero.idx])), prob)
  return(sample.mtx)
}

GeneNorm <- function(Y){
  Y = Y/rowSums(Y)*1000000
  return (log10(Y+1))
}


############### Compare Clustering Result ########################################
candidate_thres <- c(0, 0.02, 0.05, 0.1, 0.2)
datasets <-c("Usoskin","Ting","Pollen","Deng")
methods <- c("raw","rmd","sci","mag")
repeats <- 10
set.seed(2017)
for (cid in 2:length(datasets)){
	path <- paste("RealData/",datasets[cid],"_RAW.txt",sep="")
	RAW <- read.table(path,header=F)
	true_labs <- RAW[1,]
	RAW <- t(RAW[-1,])
	Counts <- round(10^RAW-1)
	N=dim(RAW)[1];P=dim(RAW)[2]
	K = max(as.integer(true_labs))
	#### Full Data ###
	res.raw <- RAW
    fwrite(data.table(t(res.raw)), paste("real_data/raw_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
  	res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  			point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("real_data/sci_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
    	res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    	fwrite(data.table(t(res.rmd)), paste("real_data/rmd_", datasets[cid], "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    } 		
	res.mag <-magic(res.raw, genes="all_genes")$result
    fwrite(data.table(t(res.mag)), paste("real_data/mag_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)	
	#### Drop Model ###
	Pro <- exp(-0.2*RAW*RAW)
	for(i in 1:repeats){
		res.raw <- RAW * matrix(rbinom(N*P, 1, as.matrix(1 - Pro)), N, P)
		res.raw <- res.raw[,colSums(res.raw)!=0]
    	fwrite(data.table(t(res.raw)), paste("real_data/drop_raw_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
  		res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  				point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  		res.sci = t(res.sci - log10(1.01))
    	fwrite(data.table(t(res.sci)), paste("real_data/drop_sci_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	for (qq in candidate_thres){
    		res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    		fwrite(data.table(t(res.rmd)), paste("real_data/drop_rmd_", datasets[cid], "_rep", i, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	} 		
		res.mag <-magic(res.raw, genes="all_genes")$result
    	fwrite(data.table(t(res.mag)), paste("real_data/drop_mag_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)	
	}

	#### Down Sampling ####
	for (i in 1:repeats){
		num.cores = 10
		registerDoParallel(cores=num.cores)
		subcounts <- Counts
		n = dim(Counts)[1]
		result <- foreach(j=1:n)%dopar%subsample(Counts,j,0.1)
		for (j in 1:n){
			subcounts[j,] = result[[j]]
		}
		res.raw <- GeneNorm(subcounts)
		res.raw <- res.raw[,colSums(res.raw)!=0]
    	fwrite(data.table(t(res.raw)), paste("real_data/down_raw_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
  		res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  				point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  		res.sci = t(res.sci - log10(1.01))
    	fwrite(data.table(t(res.sci)), paste("real_data/down_sci_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	for (qq in candidate_thres){
    		res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    		fwrite(data.table(t(res.rmd)), paste("real_data/down_rmd_", datasets[cid], "_rep", i, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	} 		
		res.mag <-magic(res.raw, genes="all_genes")$result
    	fwrite(data.table(t(res.mag)), paste("real_data/down_mag_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)	
	}
	########### Low Expression ###
	MeanExpression <- colSums(RAW)/colSums(RAW!=0)
	res.raw <- RAW[,MeanExpression<quantile(MeanExpression)[2]]
    fwrite(data.table(t(res.raw)), paste("real_data/low_raw_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
  	res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  			point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("real_data/low_sci_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
    	res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    	fwrite(data.table(t(res.rmd)), paste("real_data/low_rmd_", datasets[cid], "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    } 		
	res.mag <-magic(res.raw, genes="all_genes")$result
    fwrite(data.table(t(res.mag)), paste("real_data/low_mag_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)	
}


##################### plot the result #######################
library("ggplot2")

colors = c("#ff7f00","#377eb8","#e41a1c","#4daf4a","#984ea3")
methods <- c("raw","rmd","sci","mag")
METHODs <- c("RAW", "scRMD", "scImpute", "MAGIC")
quantiles <- c(0, 0.02, 0.05, 0.1, 0.2)
datasets <-c("Usoskin","Ting","Pollen","Deng")
repeats = 10

################## Estimated K ##################################
datasets <-c("Deng")
for (data in datasets){
	path <- paste("RealData/",data,"_RAW.txt",sep="")
	RAW <- read.table(path,header=F)
	true_labs <- RAW[1,]
	#K = max(as.integer(true_labs))
	fulldata.df <- data.frame()
	for (i in c(1:4)){
		res_tmp = data.frame(model = rep("FullData", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		if (i ==2){
			df = t(read.table(paste("real_data/", methods[i], "_", data, "_candidate0.05.txt", sep=""), header = F, sep = ","))
		}
		else{
			df = t(read.table(paste("real_data/", methods[i], "_", data, ".txt", sep=""), header = F, sep = ","))
		}
		if (data == "Usoskin"){
			tags = t(round(2^df-1))
		}
		else{
			tags = t(round(10^df-1))
		}
		sData <- CIDR(tags)
		K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		fulldata.df <- rbind(fulldata.df, res_tmp)
	}

	dropmodel.df <- data.frame()
	for (i in c(1:3)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			if (i ==2){
				df = t(read.table(paste("real_data/drop_", methods[i], "_", data, "_rep", j, "_candidate0.05.txt", sep=""), header = F, sep = ","))
			}
			else{
				df = t(read.table(paste("real_data/drop_", methods[i], "_", data, "_rep", j, ".txt", sep=""), header = F, sep = ","))
			}
			if (data == "Usoskin"){
				tags = t(round(2^df-1))
			}
			else{
				tags = t(round(10^df-1))
			}
			sData <- CIDR(tags)
			K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DropModel", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		dropmodel.df <- rbind(dropmodel.df, res_tmp)
	}

	downsampling.df <- data.frame()
	for (i in c(1:3)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			if (i ==2){
				df = t(read.table(paste("real_data/down_", methods[i], "_", data, "_rep", j, "_candidate0.05.txt", sep=""), header = F, sep = ","))
			}
			else{
				df = t(read.table(paste("real_data/down_", methods[i], "_", data, "_rep", j, ".txt", sep=""), header = F, sep = ","))
			}
			tags = t(round(10^df-1))
			sData <- CIDR(tags)
			K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DownSampling", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		downsampling.df <- rbind(downsampling.df, res_tmp)
	}

	lowexpression.df <- data.frame()
	for (i in c(1:3)){
		res_tmp = data.frame(model = rep("LowExpression", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		if (i ==2){
			df = t(read.table(paste("real_data/low_", methods[i], "_", data, "_candidate0.05.txt", sep=""), header = F, sep = ","))
		}
		else{
			df = t(read.table(paste("real_data/low_", methods[i], "_", data, ".txt", sep=""), header = F, sep = ","))
		}
		if (data == "Usoskin"){
			tags = t(round(2^df))
		}
		else{
			tags = t(round(10^df))
		}
		index = which(apply(df, 1, sum)==0)
		df[index, ] = rnorm(length(index) * dim(df)[2], 0.1, 0.1)
		sData <- CIDR(tags)
		K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		lowexpression.df <- rbind(lowexpression.df, res_tmp)
	}
	all_df = rbind(fulldata.df, dropmodel.df, downsampling.df, lowexpression.df)
	write.table(all_df, file = paste("real_data/", data, "_res.txt", sep=""), quote = F, col.names = T, row.names = F)	
}

########################################## plot result ##############################


############################ sensitivity analysis ####################################
QUANRILES = c("cutoff_0", "cutoff_0.02", "cutoff_0.05", "cutoff_0.1", "cutoff_0.2")
datasets <-c("Ting","Pollen","Deng")
for (data in datasets){
	path <- paste("RealData/",data,"_RAW.txt",sep="")
	RAW <- read.table(path,header=F)
	true_labs <- RAW[1,]
	K = max(as.integer(true_labs))
	fulldata.df <- data.frame()
	for (i in c(1:5)){
		res_tmp = data.frame(model = rep("FullData", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		df = t(read.table(paste("real_data/rmd_", data, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
		if (data == "Usoskin"){
			tags = t(round(2^df-1))
		}
		else{
			tags = t(round(10^df-1))
		}
		sData <- CIDR(tags)
		#K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		fulldata.df <- rbind(fulldata.df, res_tmp)
	}

	dropmodel.df <- data.frame()
	for (i in c(1:5)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			df = t(read.table(paste("real_data/drop_rmd_", data, "_rep", j, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
			if (data == "Usoskin"){
				tags = t(round(2^df-1))
			}
			else{
				tags = t(round(10^df-1))
			}
			sData <- CIDR(tags)
			#K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DropModel", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		dropmodel.df <- rbind(dropmodel.df, res_tmp)
	}

	downsampling.df <- data.frame()
	for (i in c(1:5)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			df = t(read.table(paste("real_data/down_rmd_", data, "_rep", j, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
			tags = t(round(10^df-1))
			sData <- CIDR(tags)
			#K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DownSampling", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		downsampling.df <- rbind(downsampling.df, res_tmp)
	}

	lowexpression.df <- data.frame()
	for (i in c(1:5)){
		res_tmp = data.frame(model = rep("LowExpression", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		df = t(read.table(paste("real_data/low_rmd_", data, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
		if (data == "Usoskin"){
			tags = t(round(2^df))
		}
		else{
			tags = t(round(10^df))
		}
		index = which(apply(df, 1, sum)==0)
		df[index, ] = rnorm(length(index) * dim(df)[2], 0.1, 0.1)
		sData <- CIDR(tags)
		#K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		lowexpression.df <- rbind(lowexpression.df, res_tmp)
	}
	all_df = rbind(fulldata.df, dropmodel.df, downsampling.df, lowexpression.df)
	write.table(all_df, file = paste("real_data/", data, "_sensitivity.txt", sep=""), quote = F, col.names = T, row.names = F)	
}

######
library(ggplot2)
colors = c("#ff7f00","#377eb8","#e41a1c","#4daf4a","#984ea3")
datasets <-c("Usoskin","Ting","Pollen","Deng")

for (data in datasets){
	df <- read.table(paste(data, "_sensitivity.txt", sep=""), header = T)
	df$xx = paste(df$model, df$clusm)
	pdf(paste(data,"_sensitivity.pdf",sep=""),height=5,width=16)
		p<-ggplot(data=df, aes(x=xx, y=ari, fill=method)) +
			geom_bar(stat="identity", position=position_dodge())+ 
			scale_fill_manual(values=colors)+
			ggtitle(data)+theme(plot.title = element_text(hjust = 0.5))+
  			theme_minimal()+geom_errorbar(aes(ymin=ari-arisd, ymax=ari+arisd), width=.2,position=position_dodge(.9)) 		
  		print(p)

	dev.off()
}

############################# Paul dataset ###############################################
paul_raw <- read.table("RealData/paul_raw.txt")


LibSize <- colSums(paul_raw)
paul_filter <- t(paul_raw)/LibSize * median(LibSize)
index = which(colSums(paul_filter > 0)/dim(paul_filter)[1] >= 0.05)
paul_filter <- paul_filter[, index]

res.raw <- GeneNorm(paul_filter)
res.rmd_ori <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw > 0], 0.05))
res.rmd.L <- data.frame(res.rmd_ori$L)
colnames(res.rmd.L) = colnames(res.raw)
res.rmd.expr <- data.frame(res.rmd_ori$exprs)

facs <- read.table("RealData/paul_info.txt", header = T)

paul_filter_p <- cbind(paul_filter, facs[, 12:13])
res.raw_p <- cbind(res.raw, facs[, 12:13])
res.rmd.L_p <- cbind(res.rmd.L, facs[, 12:13])
res.rmd.expr_p <- cbind(res.rmd.expr, facs[, 12:13]
index = intersect(which(facs$CD34_measurement >= 1), which(facs$FcgR3_measurement >= 1))

facs.plot <- ggplot(data.frame(facs[index, ]), aes(y=log(FcgR3_measurement), x=log(CD34_measurement))) + 
			geom_point(color = "blue") +
  			geom_smooth(method=lm, se=FALSE, color="red") +
  			ggtitle("FACS") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Cd34") + ylab("Fcgr3") +
  			theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 20)) +
  			annotate(geom="text", x=2.5, y=8, label=paste("r = ", round(cor(log(facs[index, ]$FcgR3_measurement), log(facs[index, ]$CD34_measurement)), 3)), color="black", size = 5)

raw.plot <- ggplot(data.frame(res.raw_p[index, ]), aes(y=Fcgr3, x=Cd34)) + 
			geom_point(color = "blue") +
  			geom_smooth(method=lm, se=FALSE, color="red") +
  			ggtitle("mRNA") + theme(plot.title = element_text(hjust = 0.5)) +
  			theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 20)) +
  			annotate(geom="text", x=2.5, y=4, label=paste("r = ", round(cor(res.raw_p[index, ]$Cd34, res.raw_p[index, ]$Fcgr3), 3)), color="black", size = 5)

rmd.L.plot <- ggplot(data.frame(res.rmd.L_p[index, ]), aes(y=Fcgr3, x=Cd34)) + 
			geom_point(color = "blue") +
  			geom_smooth(method=lm, se=FALSE, color="red") +
  			ggtitle("scRMD L") + theme(plot.title = element_text(hjust = 0.5)) +
  			theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 20)) +
  			annotate(geom="text", x=2, y=4, label=paste("r = ", round(cor(res.rmd.L_p[index, ]$Cd34, res.rmd.L_p[index, ]$Fcgr3), 3)), color="black", size = 5)

rmd.expr.plot <- ggplot(data.frame(res.rmd.expr_p[index, ]), aes(y=Fcgr3, x=Cd34)) + 
			geom_point(color = "blue") +
  			geom_smooth(method=lm, se=FALSE, color="red") +
  			ggtitle("scRMD exprs") + theme(plot.title = element_text(hjust = 0.5)) +
  			theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 20)) +
  			annotate(geom="text", x=2.5, y=4, label=paste("r = ", round(cor(res.rmd.expr_p[index, ]$Cd34, res.rmd.expr_p[index, ]$Fcgr3), 3)), color="black", size = 5)

pdf("paul_correlation.pdf",width=16,height=4)
multiplot(facs.plot, raw.plot, rmd.L.plot, rmd.expr.plot, cols=4)
dev.off()

########################################################### Zeigenhain #########################################################################

Methods <- c("CELseq","DropSeq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")
threshold = 0.05
Zeig <- read.table("RealData/Zeig_complete.txt",header=T)
ERCC <- read.table("RealData/ERCC.txt",header=F)
for (i in 1:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	else{	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	
	correlation <- data.frame(RAW=numeric(length(ind)),scRMD=numeric(length(ind)),scImpute=numeric(length(ind)),MAGIC=numeric(length(ind)))
	data <- Zeig[,ind]
	data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	logdata <- log10(data+1)
	res.rmd <- rmd(as.matrix(logdata))
	plist = get_mix_parameters(logdata + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = logdata+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	res.sci = res.sci - log10(1.01)
	res.mag <- t(magic(t(as.matrix(logdata))))

	ercc_index = grep("ERCC",rownames(logdata))
	index <- GetIndex(rownames(logdata)[ercc_index],as.character(ERCC[,2]))
	for (j in 1:length(ind)){
		correlation[j,1]=cor(logdata[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,2]=cor(res.rmd$exprs[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,3]=cor(res.sci[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,4]=cor(res.mag[ercc_index,j],log10(ERCC[index,4]))
	}
	corre <- melt(correlation)
	colnames(corre)=c("method","correlation")
	pdf(paste("barplot",Methods[i],".pdf",sep=""),height=4,width=4)
	p <- ggplot(corre, aes(x=method, y=correlation,fill=method))+ 
  	geom_boxplot(outlier.colour="red", outlier.shape=8,
                outlier.size=4) +scale_fill_manual(values=c("#377eb8","#e41a1c","#4daf4a","#984ea3"))
	print(p)
	dev.off()

}


##### Wilcoxon signed-rank test ####

Wilcoxon <- data.frame(RAW.scRMD=numeric(5),RAW.scImpute=numeric(5),
	RAW.MAGIC=numeric(5),scRMD.scImpute=numeric(5),
	scRMD.MAGIC=numeric(5),scImpute.MAGIC=numeric(5))

Methods <- c("CELseq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")
threshold = 0.05
Zeig <- read.table("RealData/Zeig_complete.txt",header=T)
ERCC <- read.table("RealData/ERCC.txt",header=F)
for (i in 1:5){
	if (i==5){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	if (i<5){	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	
	correlation <- data.frame(RAW=numeric(length(ind)),scRMD=numeric(length(ind)),scImpute=numeric(length(ind)),MAGIC=numeric(length(ind)))
	data <- Zeig[,ind]
	data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	logdata <- log10(data+1)
	res.rmd <- rmd(as.matrix(logdata))
	plist = get_mix_parameters(logdata + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = logdata+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	res.sci = res.sci - log10(1.01)
	res.mag <- t(magic(t(as.matrix(logdata))))

	ercc_index = grep("ERCC",rownames(logdata))
	index <- GetIndex(rownames(logdata)[ercc_index],as.character(ERCC[,2]))
	for (j in 1:length(ind)){
		correlation[j,1]=cor(logdata[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,2]=cor(res.rmd$exprs[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,3]=cor(res.sci[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,4]=cor(res.mag[ercc_index,j],log10(ERCC[index,4]))
	}
	Wilcoxon[i,1] = wilcox.test(correlation$RAW,correlation$scRMD, paired=FALSE)$p.value
	Wilcoxon[i,2] = wilcox.test(correlation$RAW,correlation$scImpute, paired=FALSE)$p.value
	Wilcoxon[i,3] = wilcox.test(correlation$RAW,correlation$MAGIC, paired=FALSE)$p.value
	Wilcoxon[i,4] = wilcox.test(correlation$scRMD,correlation$scImpute, paired=FALSE)$p.value
	Wilcoxon[i,5] = wilcox.test(correlation$scRMD,correlation$MAGIC, paired=FALSE)$p.value
	Wilcoxon[i,6] = wilcox.test(correlation$scImpute,correlation$MAGIC, paired=FALSE)$p.value
}

write.table(t(Wilcoxon),file="Wilcoxon.txt",quote=F,col.names=F)
############ Mean Error ################

Methods <- c("CELseq","DropSeq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")

threshold = 0.05
Zeig <- read.table("RealData/Zeig_complete.txt",header=T)
ERCC <- read.table("RealData/ERCC.txt",header=F)
MeanError <- data.frame(scRMD=numeric(6),scImpute=numeric(6),MAGIC=numeric(6))
for (i in 1:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	if(i!=6){	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	data <- Zeig[,ind]
	data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	logdata <- log10(data+1)
	Pro <- exp(-0.2*logdata*logdata)
	N = dim(logdata)[1]; P = dim(logdata)[2]
	set.seed(2017)
	raw <- logdata * matrix(rbinom(N*P, 1, 1 - Pro), N, P)
	drop <- which(raw==0)
	res.rmd <- rmd(as.matrix(raw))
	plist = get_mix_parameters(raw + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = raw+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	res.sci = res.sci - log10(1.01)
	res.mag <- t(magic(t(as.matrix(raw))))
	MeanError[i,1] = mean(abs((res.rmd$exprs-logdata)[drop]))
	MeanError[i,2] = mean(abs((res.sci-logdata)[drop]))
	MeanError[i,3] = mean(abs((res.mag-logdata)[drop]))
}


######### Cell Correlation ################
Methods <- c("CELseq","DropSeq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")
alpha = 0.1
ind <- grep(Methods[1],colnames(Zeig))
data <- Zeig[,ind]
N = dim(data)[2]
RowSum <- apply(data>0,1,sum)
index = which(RowSum/N>=alpha)

for (i in 2:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	if(i!=6){	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	data <- Zeig[,ind]
	N = dim(data)[2]
	RowSum <- apply(data>0,1,sum)
	index <- intersect(index,which(RowSum/N>=alpha))
}

HighZeig <- Zeig[index,]
ImputedData <- list(raw=list(),rmd=list(),sci=list(),mag=list())

for (i in 1:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(HighZeig)),grep(Methods[i-1],colnames(HighZeig)))
	}
	if(i!=6){	
		ind <- grep(Methods[i],colnames(HighZeig))
	}
	data <- HighZeig[,ind]
	#data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	raw <- log10(data+1)
	ImputedData$raw[[i]] <- raw
	ImputedData$rmd[[i]] <- rmd(as.matrix(raw))$exprs
	plist = get_mix_parameters(raw + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = raw+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	ImputedData$sci[[i]] = res.sci - log10(1.01)
	ImputedData$mag[[i]] <- t(magic(t(as.matrix(raw))))
}

correlation <- data.frame(RAW=numeric(length(15)),scRMD=numeric(length(15)),scImpute=numeric(length(15)),MAGIC=numeric(length(15))) 

k=1
for (i in 1:5){
	for (j in (i+1):6){
		correlation[k,1] = median(cor(ImputedData$raw[[i]],ImputedData$raw[[j]]))
		correlation[k,2] = median(cor(ImputedData$rmd[[i]],ImputedData$rmd[[j]]))
		correlation[k,3] = median(cor(ImputedData$sci[[i]],ImputedData$sci[[j]]))
		correlation[k,4] = median(cor(ImputedData$mag[[i]],ImputedData$mag[[j]]))
		k=k+1
	}
}


corre <- melt(correlation)
colnames(corre)=c("method","correlation")
pdf("barplot.pdf",height=4,width=3)
p <- ggplot(corre, aes(x=method, y=correlation)) + 
geom_boxplot(outlier.colour="red", outlier.shape=8,
    outlier.size=4)+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)
p
dev.off()
