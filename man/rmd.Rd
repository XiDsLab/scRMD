\name{rmd}
\alias{rmd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Roubst matrix decomposition
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rmd(Y, tau = NULL, lambda = NULL, initL = NULL, initS = NULL, initLambda = NULL, maxiter = 100, abstol = 0.001, reltol = 0.001, rho = 1, overrelax = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{initL}{
%%     ~~Describe \code{initL} here~~
}
  \item{initS}{
%%     ~~Describe \code{initS} here~~
}
  \item{initLambda}{
%%     ~~Describe \code{initLambda} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{abstol}{
%%     ~~Describe \code{abstol} here~~
}
  \item{reltol}{
%%     ~~Describe \code{reltol} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{overrelax}{
%%     ~~Describe \code{overrelax} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, tau = NULL, lambda = NULL, initL = NULL, initS = NULL, 
    initLambda = NULL, maxiter = 100, abstol = 0.001, reltol = 0.001, 
    rho = 1, overrelax = 1.5) 
{
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    Omega <- (Y != 0)
    if (is.null(lambda)) 
        lambda <- max(sqrt(n), sqrt(p)) * sd(Y[Omega])
    L.d <- svd(Y, 0, 0)$d
    econ <- ifelse(min(L.d) < lambda, 1, 0)
    if (is.null(initL)) {
        initL <- svt(Y, lambda, econ)$A.svt
    }
    if (is.null(initS)) {
        initS <- matrix(0, n, p)
        initS[!Omega] <- initL[!Omega]
    }
    if (is.null(tau)) 
        tau <- sd(Y[Omega])
    if (is.null(initLambda)) {
        initLambda <- matrix(0, n, p)
    }
    L <- initL
    S <- initS
    Z <- L
    Lambda <- initLambda
    alpha <- overrelax
    history <- list(rho = c(), s_norm = c(), r_norm = c(), tol_pri = c(), 
        tol_dual = c())
    for (k in 1:maxiter) {
        S[Omega] <- 0
        tmp <- (Y[Omega] + rho * L[Omega] - Lambda[Omega])/(1 + 
            rho)
        Z[Omega] = pmax(tmp, 0)
        index <- L[!Omega] < (1 + rho) * tau/rho + Lambda[!Omega]/rho
        tmp1 <- pmax((rho * L[!Omega] - Lambda[!Omega])/(1 + 
            rho), 0)
        tmp2 <- pmax(L[!Omega] - Lambda[!Omega]/rho - tau/rho, 
            0)
        tmpS <- pmax(tmp2 - tau, 0)
        tmpS[index] <- 0
        tmpZ <- tmp2
        tmpZ[index] = tmp1[index]
        S[!Omega] = tmpS
        Z[!Omega] = tmpZ
        Z_hat <- alpha * Z + (1 - alpha) * L
        L_old <- L
        tmp <- Z_hat + Lambda/rho
        svts <- svt(tmp, lambda/rho, econ)
        L <- svts$A.svt
        r <- svts$r
        Lambda <- Lambda + rho * (Z_hat - L)
        history$rho[k] <- rho
        history$r_norm[k] <- norm(Z - L, "F")
        history$s_norm[k] <- norm((L - L_old), "F") * rho
        history$tol_pri[k] <- sqrt(n * p) * abstol + reltol * 
            max(norm(Z, "F"), norm(L, "F"))
        history$tol_dual[k] = sqrt(p * n) * abstol + reltol * 
            norm(Lambda, "F")
        if (history$r_norm[k] < history$tol_pri[k] && history$s_norm[k] < 
            history$tol_dual[k]) 
            break
        if (history$r_norm[k] > 10 * history$s_norm[k]) {
            rho <- rho * 2
        }
        else if (history$s_norm[k] > 10 * history$r_norm[k]) {
            rho <- max(rho/2, 1e-04)
        }
    }
    exprs <- Y
    exprs[!Omega] <- L[!Omega]
    exprs[exprs < 0] <- 0
    return(list(L = L, S = S, r = r, Lambda = Lambda, exprs = exprs, 
        tau = tau, lambda = lambda, history = history))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
