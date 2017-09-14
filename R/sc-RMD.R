#' @import RSpectra
#' @import corpcor


svt <- function(A, thresh, econ = 0) {
  # econ is a flag that determines whether to use fast.svd in the "corpcor" package
  # or to use svds in the "RSpectra" pacakge
  if (econ == 0){
    res <- fast.svd(A)
    A.svt <- res$u %*% diag(pmax(res$d - thresh, 0)) %*% as.matrix(t(res$v))
  }
  else if (econ == 1) {
    A.d <- svd(A, 0, 0)$d
    k <- sum(A.d > thresh)
    res <- svds(A, k)
    if (k > 1) {
      d <- diag(pmax(res$d - thresh, 0))
    } else {
      d <- res$d
    }
    A.svt <- res$u %*% d %*% t(res$v)
  }
  return(list(A.svt=A.svt, r=k))
}

rmd <- function(Y, tau = NULL, lambda = NULL, initL = NULL, initS = NULL, initLambda = NULL, maxiter = 100,
                abstol = 1e-3, reltol = 1e-3, rho = 1, overrelax = 1.5) {
  # minimize 1/2||Y - (Z - S)||_F^2 + lambda*||L||_* + tau * ||S||_1
  # suject to Z = L, Z >= 0, P_Omega(S) = 0, P_{Omega^c}(S) >= 0
  
  # initialization
  n <- dim(Y)[1] # n and p is exchangeble actually
  p <- dim(Y)[2]
  Omega <- (Y != 0)
  if (is.null(lambda)) lambda <-  max(sqrt(n),sqrt(p))*sd(Y[Omega])
  # determine whether to use svds
  L.d <- svd(Y, 0, 0)$d
  econ <- ifelse(min(L.d) < lambda, 1, 0)
  if (is.null(initL)) {
    initL <- svt(Y, lambda, econ)$A.svt
  }
  if (is.null(initS)) {
    initS <- matrix(0, n, p)
    initS[!Omega] <- initL[!Omega]
  }
  #lambda <- lambda*sd(Y[Omega] - initL[Omega])
  if (is.null(tau)) tau <- sd(Y[Omega])
  if (is.null(initLambda)) {
    initLambda <- matrix(0, n, p)
  }
  L <- initL; S <- initS; Z <- L; Lambda <- initLambda; alpha <- overrelax
  
  # solve
  history <- list(rho = c(), s_norm = c(), r_norm = c(), tol_pri = c(), tol_dual = c())
  for (k in 1:maxiter) {
    # update Z,S
    # on Omega
    S[Omega] <- 0 
    tmp <- (Y[Omega] + rho * L[Omega] - Lambda[Omega]) / (1 + rho)
    Z[Omega] = pmax(tmp, 0)
    # on !Omega
    index <- L[!Omega] < (1+rho)*tau/rho + Lambda[!Omega] / rho
    tmp1 <- pmax((rho * L[!Omega] - Lambda[!Omega]) / (1 + rho), 0) # s = 0
    tmp2 <- pmax(L[!Omega] - Lambda[!Omega] / rho - tau / rho, 0) # s = z - tau
    tmpS <- pmax(tmp2 - tau, 0); tmpS[index] <- 0
    tmpZ <- tmp2; tmpZ[index] = tmp1[index]
    S[!Omega] = tmpS
    Z[!Omega] = tmpZ
    # overrelaxation
    Z_hat <- alpha * Z + (1 - alpha) * L;
    
    # update L
    L_old <- L;
    tmp <- Z_hat + Lambda / rho;
    svts <- svt(tmp, lambda / rho, econ)
    L <- svts$A.svt
    r <- svts$r
    
    # update Lambda
    Lambda <- Lambda + rho * (Z_hat - L)
    
    # diagnostics
    history$rho[k] <- rho;
    
    history$r_norm[k] <- norm(Z - L, "F")
    history$s_norm[k] <- norm((L - L_old), "F")*rho
    
    history$tol_pri[k] <- sqrt(n*p)*abstol + reltol*max(norm(Z, "F"), norm(L, "F"))
    history$tol_dual[k]= sqrt(p*n)*abstol + reltol*norm(Lambda,"F");
    
    if (history$r_norm[k] < history$tol_pri[k] && history$s_norm[k] < history$tol_dual[k])  break
    if (history$r_norm[k] > 10*history$s_norm[k]) {
      rho <- rho*2
    } else if (history$s_norm[k] > 10*history$r_norm[k]) {
      rho <- max(rho/2,1e-4)
    }
  }
  exprs <- Y
  exprs[!Omega] <- L[!Omega]
  exprs[exprs<0] <- 0
  return(list(L=L, S=S, r=r, Lambda = Lambda, exprs = exprs, tau = tau, lambda = lambda, history = history))
}

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

MatrixFusion <- function(true_lab,pre_lab){
  # true_lab is the true community label
  # pre_lab is the predicted community label
  r = max(true_lab)
  s = max(pre_lab)
  A = matrix(0,r,s)
  for (i in 1:r){
    for(j in 1:s){
      A[i,j] = length(intersect(which(true_lab==i),which(pre_lab==j)))
    }
  }
  return(A)
}


CalculateARI <- function(true_lab,pre_lab){
  # true_lab is the true community label
  # pre_lab is the predicted community label
  # return the adjusted rand index
  A = MatrixFusion(true_lab,pre_lab)
  a = apply(A,1,sum)
  b = apply(A,2,sum)
  n = length(true_lab)
  c1 = sum(a*(a-1))/2
  c2 = sum(b*(b-1))/2
  c0 = sum(A*(A-1))/2
  ari = (n*(n-1)*c0/2-c1*c2)/(n*(n-1)*(c1+c2)/4-c1*c2)
  return(ari)
}

CalculateNMI <- function(true_lab,pre_lab){
  # true_lab is the true community label
  # pre_lab is the predicted community label
  # return the normalized mutual information
  A = MatrixFusion(true_lab,pre_lab)
  n = sum(A)
  A = A/n
  p1 = apply(A,1,sum)
  p2 = apply(A,2,sum)
  I = 0
  H1 = 0
  H2 = 0
  for (i in 1:length(p1)){
    for (j in 1:length(p2)){
      if (A[i,j]!=0){
        I = I +A[i,j]*log(A[i,j]/p1[i]/p2[j])
      }
    }
  }

  for (i in 1:length(p1)){
    H1 = H1-p1[i]*log(p1[i])
  }
  for (i in 1:length(p2)){
    H2 = H2-p2[i]*log(p2[i])
  }

  return (2*I/(H1+H2))
}