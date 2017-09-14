\name{sSimulator}
\alias{sSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate scRNA data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sSimulator(K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, dr, sigmahomo, sigmahetero, type = "cluster")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{Kn}{
%%     ~~Describe \code{Kn} here~~
}
  \item{Ndiff}{
%%     ~~Describe \code{Ndiff} here~~
}
  \item{Nsame}{
%%     ~~Describe \code{Nsame} here~~
}
  \item{logMean}{
%%     ~~Describe \code{logMean} here~~
}
  \item{logSd}{
%%     ~~Describe \code{logSd} here~~
}
  \item{ZeroRate}{
%%     ~~Describe \code{ZeroRate} here~~
}
  \item{drbase}{
%%     ~~Describe \code{drbase} here~~
}
  \item{dr}{
%%     ~~Describe \code{dr} here~~
}
  \item{sigmahomo}{
%%     ~~Describe \code{sigmahomo} here~~
}
  \item{sigmahetero}{
%%     ~~Describe \code{sigmahetero} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (K, Kn, Ndiff, Nsame, logMean, logSd, ZeroRate, drbase, 
    dr, sigmahomo, sigmahetero, type = "cluster") 
{
    K <- K
    N <- K * Kn
    P <- Ndiff + Nsame
    Label <- rep(1:K, each = Kn)
    Z <- matrix(0, N, K)
    for (i in 1:K) {
        Z[((i - 1) * Kn + 1):(i * Kn), i] = 1
    }
    Esame <- rnorm(Nsame, logMean, logSd)
    Esame <- Esame * rbinom(Nsame, 1, 1 - ZeroRate)
    Esame = rep(Esame, K)
    dim(Esame) <- c(Nsame, K)
    if (type == "cluster") {
        Ediff <- matrix(rnorm(K * Ndiff, logMean, logSd), nrow = K, 
            Ndiff)
    }
    else if (type == "DE") {
        Ediff <- rnorm(Ndiff, logMean, logSd)
        Ediff <- matrix(Ediff, K, Ndiff, byrow = T)
        Ndiff.u <- round(Ndiff/K)
        for (k in 1:K) {
            Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] <- Ediff[k, 
                ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] * 2
        }
    }
    Ediff <- Ediff * matrix(rbinom(K * Ndiff, 1, 1 - ZeroRate), 
        K, Ndiff)
    Te <- Z \%*\% cbind(Ediff, t(Esame))
    Te[Te > 0] <- Te[Te > 0] + rnorm(sum(Te > 0), 0, Te[Te > 
        0] * sigmahetero + sigmahomo)
    Te[Te <= 0] <- rnorm(sum(Te <= 0), 0, 1) - 1.5
    Te[Te <= 0] <- 0
    Pro <- drbase * exp(-dr * Te * Te)
    De <- Te * matrix(rbinom(N * P, 1, 1 - Pro), N, P)
    sData = list(K = K, N = N, P = P, label = Label, te = Te, 
        de = De)
    return(sData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
