\name{MatrixFusion}
\alias{MatrixFusion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MatrixFusion(true_lab, pre_lab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{true_lab}{
%%     ~~Describe \code{true_lab} here~~
}
  \item{pre_lab}{
%%     ~~Describe \code{pre_lab} here~~
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
function (true_lab, pre_lab) 
{
    r = max(true_lab)
    s = max(pre_lab)
    A = matrix(0, r, s)
    for (i in 1:r) {
        for (j in 1:s) {
            A[i, j] = length(intersect(which(true_lab == i), 
                which(pre_lab == j)))
        }
    }
    return(A)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
