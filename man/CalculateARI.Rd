\name{CalculateARI}
\alias{CalculateARI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalculateARI(true_lab, pre_lab)
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
    a = apply(A, 1, sum)
    b = apply(A, 2, sum)
    n = length(true_lab)
    c1 = sum(a * (a - 1))/2
    c2 = sum(b * (b - 1))/2
    c0 = sum(A * (A - 1))/2
    ari = (n * (n - 1) * c0/2 - c1 * c2)/(n * (n - 1) * (c1 + 
        c2)/4 - c1 * c2)
    return(ari)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
