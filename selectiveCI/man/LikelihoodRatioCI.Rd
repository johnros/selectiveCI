\name{LikelihoodRatioCI}
\alias{LikelihoodRatioCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zhong and Prentice (2008) bounds for conditional CI
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LikelihoodRatioCI(x, cutoff, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
[1] Zhong, Hua, and Ross L. Prentice. 2008. “Bias-reduced Estimators and Confidence Intervals for Odds Ratios in Genome-wide Association Studies.” Biostatistics 9 (4) (October 1): 621–634. doi:10.1093/biostatistics/kxn001.
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
function (x, cutoff, alpha) 
{
    log.cond.lik <- function(theta) log(dnorm(x - theta)/(1 - 
        pnorm(cutoff - theta) + 1 - pnorm(cutoff + theta)))
    g <- function(theta) (2 * pi)^(-1/2) * (x - theta) * exp(-(x - 
        theta)^2/2) * (1 - pnorm(cutoff - theta) + 1 - pnorm(cutoff + 
        theta)) - dnorm(x - theta) * (dnorm(cutoff - theta) - 
        dnorm(cutoff + theta))
    root <- uniroot(g, c(min(x, 0), max(x, 0)))$root
    h <- function(theta) log.cond.lik(theta) - (log.cond.lik(root) - 
        qgamma(1 - alpha, shape = 1/2, scale = 2)/2)
    t1 <- root
    while (h(t1) > 0) t1 <- t1 - 1
    t2 <- root
    while (h(t2) > 0) t2 <- t2 + 1
    lower <- uniroot(h, c(t1, root))$root
    upper <- uniroot(h, c(root, t2))$root
    CI <- list(lower = lower, upper = upper)
    return(CI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
