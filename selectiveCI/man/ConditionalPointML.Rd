\name{ConditionalPointML}
\alias{ConditionalPointML}

\title{
Selective point estimation using maximum likelihood.
}
\description{
Returns a selective maximum likelihood point estimate for the mean of a Gaussian population.
}
\usage{
  ConditionalPointML(x, sigsq, cutoff,init.mu=NA,max.iter=1000,tol=10^-5)
}

\arguments{
  \item{x}{The observed value from the population.}
  \item{sigsq}{The (known) standard deviation of the population.}
  \item{cutoff}{The selection cutoff.}
  \item{init.mu}{Starting point for the Newton Raphson algorithm. If NA then the algorithm initializes from the observed value}
  \item{max.iter}{Maximal number of iteration for the Newton Raphson algorithm.}
  \item{tol}{The desired accuracy}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{Returns a numeric value of the estimated mean.}

\references{
Benjamini, Yoav, and Amit Meir. "Selective Correlations-the conditional estimators." arXiv preprint arXiv:1412.3242 (2014).‏
}

\author{
  Amit Meir
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{selectiveCI}}
}

\examples{
cutoff <- 1
sigsq <- 1
x <- sort(rnorm(1000,sd=sqrt(sigsq)))
x <- x[abs(x)>cutoff]

cond <- sapply(x,ConditionalPointML,sigsq,cutoff)

plot(x,cond,pch=".",col="red")
abline(a=0,b=1)
abline(v=c(-cutoff,cutoff),col="grey")
grid()
}
