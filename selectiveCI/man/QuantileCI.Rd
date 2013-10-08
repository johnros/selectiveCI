\name{QuantileCI}
\alias{QuantileCI}

\title{
  Conditional quantile based confidence interval
}

\description{
  Selective CI based on the quantiles of the conditonal distribution.
}

\usage{
QuantileCI(x, cutoff, alpha)
}


\arguments{
   \item{x}{ 
  The observed parameter estimate.
  }
  \item{cutoff}{ 
  The selection cutoff.
  }
  \item{alpha}{
  Desired confidence level.
  }
}

\value{
  Returns a list of upper and lower confidence limits.
}

\references{
[1] Weinstein, Asaf, William Fithian, and Yoav Benjamini. 2013. "Selection Adjusted Confidence Intervals With More Power to Determine the Sign." Journal of the American Statistical Association 108 (501): 165-176. 
}
\author{
  Assaf Weinstein
}


\seealso{
 \code{\link{ShortestAR}}
}

\examples{
  QuantileCI(2, 1, 0.05)
}


\keyword{univar}
\keyword{htest}
