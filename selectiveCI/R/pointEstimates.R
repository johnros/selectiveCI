ConditionalPointML <- function(x,sigsq,cutoff) {
  if(abs(x)<cutoff) {
    warning("abs(x) is less than the cutoff value!")
  }
  f <- function(mu,x,cutoff) dnorm(x,mean=mu,sd=sqrt(sigsq))/(1-pnorm(cutoff,mean=mu,sd=sqrt(sigsq))+pnorm(-cutoff,mean=mu,sd=sqrt(sigsq)))
  muhat <- suppressWarnings(optimize(f,interval=c(-6,6),x,cutoff,maximum=TRUE)[[1]])
  return(muhat)
}

##Testing
# ConditionalPointML(0,1,1)

# cutoff <- 1
# sigsq <- 1
# x <- sort(rnorm(1000,sd=sqrt(sigsq)))
# x <- x[abs(x)>cutoff]
# 
# cond <- sapply(x,ConditionalPointML,sigsq,cutoff)
# 
# plot(x,cond,pch=".",col="red")
# abline(a=0,b=1)
# abline(v=c(-cutoff,cutoff),col="grey")
# grid()
