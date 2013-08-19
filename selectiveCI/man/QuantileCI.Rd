\name{QuantileCI}
\alias{QuantileCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional quantile based acceptance region and confidence interval
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QuantileCI(x, cutoff, alpha)
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
function (x, cutoff, alpha) 
{
    f <- function(theta) 1 - pnorm(x - theta) - alpha/2 * Q(theta)
    lower <- uniroot(f, c(-cutoff - qnorm(1 - alpha/2), x - qnorm(1 - 
        alpha/2)))$root
    f <- function(theta) pnorm(theta - x) - (1 - alpha/2) * Q(theta)
    upper <- uniroot(f, c(x, x + qnorm(1 - alpha/2)))$root
    CI <- list(lower = lower, upper = upper)
    return(CI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
