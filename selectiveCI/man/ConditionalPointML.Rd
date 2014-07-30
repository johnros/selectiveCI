\name{ConditionalPointML}
\alias{ConditionalPointML}

\title{
Selective point estimation using maximum likelihood.
}
\description{
Returns a selective maximum likelihood point estimate for the mean of a Gaussian population.
}
\usage{
  ConditionalPointML(x, sigsq, cutoff)
}

\arguments{
  \item{x}{The observed value from the population.}
  \item{sigsq}{The (known) standard deviation of the population.}
  \item{cutoff}{The selection cutoff.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{Returns a numeric value of the estimated mean.}

%\references{
%% ~put references to the literature/web site here ~
%}

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
