\name{ModifiedPrattAR}
\alias{ModifiedPrattAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ModifiedPrattAR(theta, r, cutoff, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
function (theta, r, cutoff, alpha) 
{
    Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - 
        theta)
    f <- function(theta) (pnorm(cutoff + theta) - pnorm(cutoff - 
        theta)) - (1 - alpha) * Q(theta)
    theta1 <- uniroot(f, c(0, cutoff + qnorm(1 - alpha)))$root
    f <- function(theta) pnorm(cutoff + r * ShortestAR(theta, 
        cutoff, alpha)$l - theta) - pnorm(cutoff - theta) - (1 - 
        alpha) * Q(theta)
    thetatilde1 <- uniroot(f, c(0, theta1))$root
    thetatilde2 <- uniroot(f, c(theta1, r * 2 * qnorm(1 - alpha/2)))$root
    if (theta < 0) {
        tmp <- ModifiedPrattAR(-theta, r, cutoff, alpha)
        A <- -rev(tmp)
    }
    else {
        if (theta == 0) {
            ll <- -qnorm(1 - alpha/2 * Q(0))
            ul <- -cutoff
            lr <- cutoff
            ur <- qnorm(1 - alpha/2 * Q(0))
        }
        if (0 < theta && theta < thetatilde1) {
            l <- ShortestAR(theta, cutoff, alpha)$l
            f <- function(x) 1 - pnorm(theta - x) + 1 - pnorm(cutoff + 
                r * l - (-cutoff - x) - theta) - alpha * Q(theta)
            atilde1 <- uniroot(f, c(theta - qnorm(1 - alpha/2 * 
                Q(theta)), -cutoff))$root
            ll <- atilde1
            ul <- -cutoff
            lr <- cutoff
            ur <- cutoff + r * l - (-cutoff - atilde1)
        }
        if (thetatilde1 <= theta && theta < thetatilde2) {
            l <- ShortestAR(theta, cutoff, alpha)$l
            f <- function(x) pnorm(x + r * l - theta) - pnorm(x - 
                theta) - (1 - alpha) * Q(theta)
            atilde2 <- uniroot(f, c(cutoff, theta - qnorm((1 - 
                alpha) * Q(theta))))$root
            ll <- NA
            ul <- NA
            lr <- atilde2
            ur <- atilde2 + r * l
        }
        if (thetatilde2 <= theta) {
            l <- ShortestAR(theta, cutoff, alpha)$l
            f <- function(x) pnorm(x + r * l - theta) - pnorm(x - 
                theta) - (1 - alpha) * Q(theta)
            atilde2 <- uniroot(f, c(theta - r * l/2, theta - 
                qnorm((1 - alpha) * Q(theta))))$root
            ll <- NA
            ul <- NA
            lr <- atilde2
            ur <- atilde2 + r * l
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
