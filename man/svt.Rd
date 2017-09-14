\name{svt}
\alias{svt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
svt(A, thresh, econ = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{thresh}{
%%     ~~Describe \code{thresh} here~~
}
  \item{econ}{
%%     ~~Describe \code{econ} here~~
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
function (A, thresh, econ = 0) 
{
    if (econ == 0) {
        res <- fast.svd(A)
        A.svt <- res$u \%*\% diag(pmax(res$d - thresh, 0)) \%*\% 
            as.matrix(t(res$v))
    }
    else if (econ == 1) {
        A.d <- svd(A, 0, 0)$d
        k <- sum(A.d > thresh)
        res <- svds(A, k)
        if (k > 1) {
            d <- diag(pmax(res$d - thresh, 0))
        }
        else {
            d <- res$d
        }
        A.svt <- res$u \%*\% d \%*\% t(res$v)
    }
    return(list(A.svt = A.svt, r = k))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
