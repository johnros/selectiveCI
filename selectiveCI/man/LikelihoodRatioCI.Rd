\name{LikelihoodRatioCI}
\alias{LikelihoodRatioCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zhong and Prentice (2008) bounds for conditional CI
}
\description{
  Implements the biased reduced confidence intervals in the reference.
}

\usage{
  LikelihoodRatioCI(x, cutoff, alpha)
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
[1] Zhong, Hua, and Ross L. Prentice. 2008. "Bias-reduced Estimators and Confidence Intervals for Odds Ratios in Genome-wide Association Studies." Biostatistics 9 (4) (October 1): 621-634. doi:10.1093/biostatistics/kxn001.
}

\author{Assaf Weinstein}

\seealso{
 \code{\link{ShortestAR}}
}
\examples{
    LikelihoodRatioCI(x=2, cutoff=1, alpha=0.05)
}

\keyword{univar}
\keyword{htest}
