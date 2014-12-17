ConditionalPointML <- function(x,cutoff,sigsq,init.mu=NA,max.iter=1000,tol=10^-5) {
  comp.power <- function(mu,sig,cutoff) 1-pnorm(cutoff,mean=mu,sd=sig)+pnorm(-cutoff,mean=mu,sd=sig)
  phi <- function(c) dnorm(c,mean=muhat,sd=sig)
  
  comp.grad <- function(x,mu,sig,cutoff) {
    Qc <- comp.power(mu,sig,cutoff)
    (x-mu)/sig^2 - (phi(cutoff)-phi(-cutoff))/Qc
  }
  comp.hess <- function(x,mu,sig,cutoff) {
    Qc <- comp.power(mu,sig,cutoff)
    hess <- Qc*(phi(cutoff)*(cutoff-mu)/sig^2+phi(-cutoff)*(cutoff+mu)/sig^2) + 
      (phi(cutoff)-phi(-cutoff))^2
    hess <- -hess/Qc^2
    hess <- hess-1/sig^2
    return(hess)
  }
  
  ## If x is a vector then replace x with the mean and the variance with the variance of the mean
  x <- mean(x)
  sig <- sqrt(sigsq/length(x))
  
  if(abs(x)<cutoff) {
    warning("abs(x) is less than the cutoff value")
  }
  
  if(is.na(init.mu)) init.mu <- x
  muhat <- init.mu
  
  for(i in 1:max.iter) {
    old.muhat <- muhat
    
    grad <- comp.grad(x,muhat,sig,cutoff)
    hess <- comp.hess(x,muhat,sig,cutoff)
    
    muhat <- muhat - grad/hess
    
    if(abs(muhat-old.muhat)<tol) break
  }
  
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
