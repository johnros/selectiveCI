\name{QuasiConventionalAR}
\alias{QuasiConventionalAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QuasiConventionalAR(theta, lambda, cutoff, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
function (theta, lambda, cutoff, alpha) 
{
    Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - 
        theta)
    f <- function(theta) (pnorm(cutoff + theta) - pnorm(cutoff - 
        theta)) - (1 - alpha) * Q(theta)
    theta1 <- uniroot(f, c(0, cutoff + qnorm(1 - alpha)))$root
    f <- function(theta) 1 - pnorm(cutoff - theta) - (1 - alpha) * 
        Q(theta)
    thetastar <- uniroot(f, c(0, cutoff + qnorm(1 - alpha)))$root
    f <- function(theta) 1 + lambda * (1 - dnorm(cutoff + theta)/dnorm(qnorm(2 - 
        pnorm(cutoff + theta) - alpha * Q(theta))))
    b <- thetastar
    while (is.na(f(b)) == T) {
        b <- b + 1e-05
    }
    thetaprime1 <- uniroot(f, c(b, theta1))$root
    if (theta < 0) {
        tmp <- QuasiConventionalAR(-theta, lambda, cutoff, alpha)
        A <- -rev(tmp)
    }
    else {
        if (theta == 0) {
            ll <- -qnorm(1 - alpha/2 * Q(0))
            ul <- -cutoff
            lr <- cutoff
            ur <- qnorm(1 - alpha/2 * Q(0))
        }
        if (0 < theta && theta < thetaprime1) {
            f <- function(d) 1 + lambda * (1 - dnorm(cutoff + 
                d + theta)/dnorm(qnorm(2 - pnorm(cutoff + d + 
                theta) - alpha * Q(theta))))
            dunderbar <- max(-cutoff - theta + qnorm(1 - alpha * 
                Q(theta)), 0)
            dbar <- -cutoff - ShortestAR(theta, cutoff, alpha)$A[1]
            dstar <- uniroot(f, c(dunderbar + 1e-05, dbar))$root
            aprime1 <- theta - cutoff + qnorm(1 - pnorm(cutoff + 
                dstar + theta) + 1 - alpha * Q(theta))
            ll <- -cutoff - dstar
            ul <- -cutoff
            lr <- cutoff
            ur <- cutoff + aprime1
        }
        if (thetaprime1 <= theta && theta < theta1) {
            aprime2 <- theta - cutoff + qnorm(pnorm(cutoff - 
                theta) + (1 - alpha) * Q(theta))
            ll <- NA
            ul <- NA
            lr <- cutoff
            ur <- cutoff + aprime2
        }
        if (theta > theta1) {
            A <- ShortestAR(theta, cutoff, alpha)$A
            ll <- A[1]
            ul <- A[2]
            lr <- A[3]
            ur <- A[4]
        }
        A <- c(ll, ul, lr, ur)
    }
    return(A)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
