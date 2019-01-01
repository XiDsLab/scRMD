sSimulator <- function(K, Kn, Ndiff, Nsame, logMean, logSd,  
                       ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "cluster"){
  # K is #cell type, Kn is #cells in each cell type, Nsame is #genes with same expression profiles,
  # Ndiff is #DE genes. ZeroRate is used to simulate true zeros. drbase and dr determines dropout basline and slope
  # sigmahomo and sigmahetero determines noise, type determines how to simulate DE genes
  K <-  K
  N <-  K * Kn
  P <-  Ndiff + Nsame
  Label <- rep(1:K, each = Kn)
  Z <- matrix(0,N,K)
  for (i in 1:K){
    Z[((i-1)*Kn+1):(i*Kn),i] = 1
  }
  Esame <- rnorm(Nsame,logMean,logSd)
  Esame <- Esame * rbinom(Nsame, 1, 1 - ZeroRate)
  Esame = rep(Esame, K)
  dim(Esame) <- c(Nsame, K) 
  if (type == "cluster") {
    Ediff <- matrix(rnorm(K * Ndiff, logMean, logSd), nrow = K, Ndiff)
  } else if (type == "DE") {
    Ediff <- rnorm(Ndiff, logMean, logSd)
    Ediff <- matrix(Ediff, K, Ndiff, byrow = T)
    Ndiff.u <- round(Ndiff/K)
    for (k in 1:K) {
      Ediff[k,((k-1) * Ndiff.u + 1):(k*Ndiff.u)] <- Ediff[k,((k-1) * Ndiff.u + 1):(k*Ndiff.u)] * 2#runif(Ndiff.u, 1.5, 2.5)
    }
  }
  Ediff <- Ediff * matrix(rbinom(K * Ndiff, 1, 1 - ZeroRate), K, Ndiff)
  Te <- Z %*% cbind(Ediff, t(Esame))
  Te[Te > 0] <- Te[Te > 0] + rnorm(sum(Te>0), 0, Te[Te > 0] * sigmahetero + sigmahomo)
  Te[Te <= 0] <- rnorm(sum(Te <= 0), 0, 1)-1.5 # some genes are never expressed, add some noise
  Te[Te <= 0] <- 0
  Pro <- drbase * exp(-dr * Te*Te)
  De <- Te * matrix(rbinom(N*P, 1, 1 - Pro), N, P)
  sData = list(K=K,N=N,P=P,label=Label,te=Te,de=De)
  return(sData)
}