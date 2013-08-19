\name{QuasiConventionalCI}
\alias{QuasiConventionalCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QuasiConventionalCI(x, sigsq, lambda, cutoff, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sigsq}{
%%     ~~Describe \code{sigsq} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
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
function (x, sigsq, lambda, cutoff, alpha) 
{
    if (x < 0) {
        tmp <- QuasiConventionalCI(-x, sigsq, lambda, cutoff, 
            alpha)
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
        f <- function(theta) 1 - pnorm(cutoff - theta) - (1 - 
            alpha) * Q(theta)
        thetastar <- uniroot(f, c(0, cutoff + qnorm(1 - alpha)))$root
        f <- function(theta) 1 + lambda * (1 - dnorm(cutoff + 
            theta)/dnorm(qnorm(2 - pnorm(cutoff + theta) - alpha * 
            Q(theta))))
        b <- thetastar
        while (is.na(f(b)) == T) {
            b <- b + 1e-05
        }
        thetaprime1 <- uniroot(f, c(b, theta1))$root
        dmin <- function(theta) max(-cutoff - theta + qnorm(1 - 
            alpha * Q(theta)), 0)
        dmax <- function(theta) -cutoff - ShortestAR(theta, cutoff, 
            alpha)$A[1]
        zalphahalf <- qnorm(1 - 0.5 * alpha * Q(0))
        dminzero <- dmin(0)
        dmaxzero <- dmax(0)
        f <- function(d) 1 + lambda * (1 - dnorm(cutoff + d)/dnorm(qnorm(2 - 
            pnorm(cutoff + d) - alpha * Q(0))))
        xprime1 <- uniroot(f, c(dminzero, dmaxzero))$root + cutoff
        xprime2 <- qnorm(2 - pnorm(xprime1) - alpha * Q(0))
        xprime3 <- thetaprime1 + qnorm(pnorm(cutoff - thetaprime1) + 
            (1 - alpha) * Q(thetaprime1))
        is.neg <- 0
        if (x < 0) 
            is.neg <- 1
        x <- abs(x)
        if (cutoff < x && x < xprime1) {
            d <- x - cutoff
            f <- function(theta) 1 + lambda * (1 - dnorm(cutoff + 
                d + theta)/dnorm(qnorm(2 - pnorm(cutoff + d + 
                theta) - alpha * Q(theta))))
            g <- function(theta) dmin(theta) - (x - cutoff)
            if (x - cutoff > dmin(0)) {
                thetamin <- 0
            }
            else {
                g <- function(theta) dmin(theta) - (x - cutoff)
                thetamin <- uniroot(g, c(0, thetastar + 0.001))$root
            }
            leftend <- thetamin + 1e-05
            while (is.na(suppressWarnings(f(leftend)))) leftend <- leftend + 
                1e-05 * 5
            lower <- uniroot(f, c(leftend, thetaprime1 + 0.1))$root
            lower <- -lower
        }
        if (xprime1 <= x && x < zalphahalf) {
            lower <- 0
        }
        if (zalphahalf <= x && x < xprime2) {
            lower <- 0
        }
        if (xprime2 <= x && x < xprime3) {
            g <- function(theta) -cutoff - theta + qnorm(2 - 
                pnorm(x - theta) - alpha * Q(theta))
            f <- function(theta) 1 + lambda * (1 - dnorm(cutoff + 
                g(theta) + theta)/dnorm(qnorm(2 - pnorm(cutoff + 
                g(theta) + theta) - alpha * Q(theta))))
            h <- function(theta) 1 - pnorm(x - theta) - alpha * 
                Q(theta)
            thetamax <- uniroot(h, c(0, x))$root
            lower <- uniroot(f, c(0 - 0.1, thetamax - 0.1))$root
        }
        if (x >= xprime3) {
            lower <- ShortestCI(x, cutoff, alpha)$lower
        }
        f <- function(theta) 2 * pnorm(theta - cutoff) - 1 - 
            (1 - alpha) * Q(theta)
        theta2 <- uniroot(f, c(theta1, cutoff + qnorm(1 - alpha/2)))$root
        f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - 
            alpha) * Q(theta)
        upper <- uniroot(f, c(theta2, x + 1.1 * qnorm(1 - alpha/2)))$root
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
