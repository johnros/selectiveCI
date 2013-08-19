\name{ModifiedPrattCI}
\alias{ModifiedPrattCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ModifiedPrattCI(x, sigsq, r, cutoff, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sigsq}{
%%     ~~Describe \code{sigsq} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
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
function (x, sigsq, r, cutoff, alpha) 
{
    if (x < 0) {
        tmp <- ModifiedPrattCI(-x, sigsq, r, cutoff, alpha)
        lower <- -tmp$upper
        upper <- -tmp$lower
    }
    else {
        x <- x/sqrt(sigsq)
        cutoff <- cutoff/sqrt(sigsq)
        Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - 
            pnorm(cutoff - theta)
        f <- function(theta) (pnorm(cutoff + theta) - pnorm(cutoff - 
            theta)) - (1 - alpha) * Q(theta)
        theta1 <- uniroot(f, c(0, cutoff + qnorm(1 - alpha)))$root
        f <- function(theta) pnorm(cutoff + r * ShortestAR(theta, 
            cutoff, alpha)$l - theta) - pnorm(cutoff - theta) - 
            (1 - alpha) * Q(theta)
        thetatilde1 <- uniroot(f, c(0, theta1))$root
        thetatilde2 <- uniroot(f, c(theta1, r * 2 * qnorm(1 - 
            alpha/2)))$root
        zalphahalf <- qnorm(1 - 0.5 * alpha * Q(0))
        lzero <- ShortestAR(0, cutoff, alpha)$l
        f <- function(x) 1 - pnorm(x) + 1 - pnorm(2 * cutoff + 
            r * lzero - x) - 2 * alpha * (1 - pnorm(cutoff))
        xtilde1 <- uniroot(f, c(cutoff, zalphahalf))$root
        xtilde2 <- uniroot(f, c(zalphahalf, cutoff + r * 2 * 
            qnorm(1 - alpha/2)))$root
        ltilde1 <- ShortestAR(thetatilde1, cutoff, alpha)$l
        is.neg <- 0
        if (x < 0) 
            is.neg <- 1
        x <- abs(x)
        if (cutoff < x && x < xtilde1) {
            f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - 
                x + 2 * cutoff + r * ShortestAR(theta, cutoff, 
                alpha)$l) - alpha * Q(theta)
            lower <- uniroot(f, c(-thetatilde1 - 0.001, 0))$root
        }
        if (xtilde1 <= x && x < zalphahalf) {
            lower <- 0
        }
        if (zalphahalf <= x && x < xtilde2) {
            lower <- 0
        }
        if (xtilde2 <= x && x < cutoff + r * ltilde1) {
            f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - 
                x + 2 * cutoff + r * ShortestAR(theta, cutoff, 
                alpha)$l) - alpha * Q(theta)
            lower <- uniroot(f, c(0, thetatilde1))$root
        }
        if (x > cutoff + r * ltilde1) {
            f <- function(theta) pnorm(x - theta) - pnorm(x - 
                r * ShortestAR(theta, cutoff, alpha)$l - theta) - 
                (1 - alpha) * Q(theta)
            m <- optimize(f, c(thetatilde1, x), maximum = T)$maximum
            lower <- uniroot(f, c(thetatilde1, m))$root
        }
        f <- function(theta) pnorm(x + r * ShortestAR(theta, 
            cutoff, alpha)$l - theta) - pnorm(x - theta) - (1 - 
            alpha) * Q(theta)
        m <- optimize(f, c(x, x + r * 2 * qnorm(1 - alpha/2)), 
            maximum = T)$maximum
        upper <- uniroot(f, c(0, m))$root
        lower <- sqrt(sigsq) * lower
        upper <- sqrt(sigsq) * upper
    }
    CI <- list(lower = lower, upper = upper)
    return(CI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
