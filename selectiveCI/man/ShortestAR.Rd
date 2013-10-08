\name{ShortestAR}
\alias{ShortestAR}
\alias{ShortestCI}
\alias{QuasiConventionalAR}
\alias{QuasiConventionalCI}
\alias{ModifiedPrattAR}
\alias{ModifiedPrattCI}

\title{
Conditional acceptance region and confidence interval
}

\description{
Computes acceptance regions and confidence intervals for selected location parameters as descrigbed in the references.
}

\usage{
ShortestAR(theta, cutoff, alpha)

ShortestCI(x, sigsq, cutoff, alpha)

QuasiConventionalAR(theta, lambda, cutoff, alpha)

QuasiConventionalCI(x, sigsq, lambda ,cutoff ,alpha)

ModifiedPrattAR(theta, ratio, cutoff, alpha)

ModifiedPrattCI(x, sigsq, ratio, cutoff, alpha)

}


\arguments{
  \item{theta}{
  The location parameter of the unconditional underlying distribution.
}
  \item{cutoff}{ 
  The selection cutoff.
}
  \item{alpha}{
  The area under the density or the confidence level of the CI.
}
  \item{x}{
  Observed value of the parameter estimate X (naturaly abs(X)>c).
  }
  \item{sigsq}{
  The variance of X.
  }
}

\details{
  \code{\link{ShortestCI}},  \code{\link{ModifiedPrattCI}} and \code{\link{QuasiConventionalCI}} will return confidence region for \eqn{\theta} conditional on \eqn{|x|>cutoff}, i.e., conditional on \eqn{\theta} being selected.

  \code{\link{ShortestAR}}, \code{\link{ModifiedPrattAR}}, and \code{\link{QuasiConventionalAR}} will return intervals containing \eqn{1-\alpha} probaility from the distribution of \eqn{X s.t. |X|>c}.
  
  More details can be found in the reference.
}

\value{
Confidence interval functions will return the upper and lower bound of the confidence interval.
  
  Acceptance regions functions will return a list containing:
  \item{a}{A \eqn{1-\alpha} confidence set which is the union of two intervals.}
  \item{length}{The total length of the intervals}  
}

\references{
[1] Weinstein, Asaf, William Fithian, and Yoav Benjamini. 2013. "Selection Adjusted Confidence Intervals With More Power to Determine the Sign." Journal of the American Statistical Association 108 (501): 165-176. 
}

\author{
Jonathan Rosenblatt and Asaf Weinstein
}

\examples{
  ShortestAR(theta=1, cutoff=1, alpha=0.05)
  ShortestCI(x=1.9, sigsq=1, cutoff=1, alpha=0.01) 
}
\keyword{univar}
\keyword{htest}
