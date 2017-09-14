\name{CalculateNMI}
\alias{CalculateNMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate normalized mutual information
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalculateNMI(true_lab, pre_lab)
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
    A = MatrixFusion(true_lab, pre_lab)
    n = sum(A)
    A = A/n
    p1 = apply(A, 1, sum)
    p2 = apply(A, 2, sum)
    I = 0
    H1 = 0
    H2 = 0
    for (i in 1:length(p1)) {
        for (j in 1:length(p2)) {
            if (A[i, j] != 0) {
                I = I + A[i, j] * log(A[i, j]/p1[i]/p2[j])
            }
        }
    }
    for (i in 1:length(p1)) {
        H1 = H1 - p1[i] * log(p1[i])
    }
    for (i in 1:length(p2)) {
        H2 = H2 - p2[i] * log(p2[i])
    }
    return(2 * I/(H1 + H2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
