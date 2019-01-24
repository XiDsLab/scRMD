#'imputation by robust matrix decomposition
#'
#' @param Y A single cell RNA data matrix;
#' rows representing cells.
#' @param tau Tuning parameter to penalize the sparsity of S;
#' @param lambda Tuning parameter to penalize the row rank of L;
#' @param initL The initionlization of L;
#' @param initS The initionlization of S;
#' @param initLambda The initionlization of Lambda;
#' @param maxiter maxmium iteration of algorithm;
#' @param candidate the cutoff for candidate drop out;
#' @export
#' @import RSpectra
#' @import corpcor
#' @author Chong Chen, \email{cheung1990@126.com}

rmd <- function(Y, tau = NULL, lambda = NULL, initL = NULL, initS = NULL, initLambda = NULL, maxiter = 100,
                abstol = 1e-3, reltol = 1e-3, rho = 1, overrelax = 1.5, candidate = 0.05, econ = 1) {
  # minimize 1/2||Y - (Z - S)||_F^2 + lambda*||L||_* + tau * ||S||_1
  # suject to Z = L, Z >= 0, P_Omega(S) = 0, P_{Omega^c}(S) >= 0
  
  # initialization
  n <- dim(Y)[1] # n and p is exchangeble actually
  p <- dim(Y)[2]
  if (min(p, n) > 1000){
    econ = 0
  }
  Omega <- (Y > candidate)
  if (is.null(lambda)) {
    lambda <-  (sqrt(n) + sqrt(p)) * sd(Y)
  }
  # determine whether to use svds
  # L.d <- svd(Y, 0, 0)$d
  #econ <- ifelse(min(L.d) < lambda, 1, 0)
  if (is.null(initL)) {
    initL <- svt(Y, lambda, econ)$A.svt
  }
  if (is.null(initS)) {
    initS <- matrix(0, n, p)
    initS[!Omega] <- initL[!Omega]
  }
  #lambda <- lambda*sd(Y - initL)
  if (is.null(tau)) {
    tau <- sd(Y[Omega] - initL[Omega])
  }
  if (is.null(initLambda)) {
    initLambda <- matrix(0, n, p)
  }
  L <- initL; S <- initS; Z <- L; Lambda <- initLambda; alpha <- overrelax
  
  # solve
  history <- list(rho = c(), s_norm = c(), r_norm = c(), tol_pri = c(), tol_dual = c())
  for (k in 1:maxiter) {
  print(k)
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
    L[L<0] = 0
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
  exprs[S>0] <- L[S>0]
  exprs[exprs<0] <- 0
  return(list(L=L, S=S, r=r, Lambda = Lambda, exprs = exprs, tau = tau, lambda = lambda, history = history))
}