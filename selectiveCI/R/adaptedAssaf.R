#CCI functions
#commentary added Aug 16th, 2013.

# The following functions implement the three 
# methods suggested in the paper to construct conditional 
# CIs for theta. The first method is "Shortest-AR", 
# the second is "Modified Pratt" and 
# the third is "Quasi-Conventional". 
# There are two functions corresponding to each of the 
# three methods: The first computes the acceptance 
# regions (AR), and the second computes the CI 
# (that are obtained when inverting the ARs and 
# taking the convex hull).


# inputs are: theta (value of location parameter), cutoff (positive threshold) and alpha (level of the test)
# Note: ARs are constructed for sigsq=1; This is enough to obtain a CI for a general sigsq (proper modification is made in the CI function)
ShortestAR <- function(theta,cutoff,alpha){
  Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)	
  #compute useful quantities
  f <- function(theta) ( pnorm(cutoff + theta) - pnorm(cutoff - theta) ) - (1 - alpha) * Q(theta)
  theta1 <- uniroot(f,c(0,cutoff + qnorm(1 - alpha))) $ root
  f <- function(theta) 2*pnorm(theta - cutoff) -1 - (1 - alpha)* Q(theta)
  theta2 <- uniroot( f,c( theta1,cutoff + qnorm(1 - alpha/2) ) ) $ root      
  # theta2 - cutoff cannot be greater than qnorm(1 - alpha/2)!
  
  if (theta<0) {
    tmp <- Recall(theta=-theta, cutoff=cutoff, alpha=alpha)
    A <- -rev(tmp$A)
    l <- tmp$l
  }
  else{
    #compute ends of the AR
    if (theta == 0) {
      ll <- -qnorm( 1 - alpha/2 * Q(0) )
      ul <- -cutoff
      lr <- cutoff
      ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
    
    if (0 < theta && theta < theta1) {
      ll <- theta - qnorm( 1 - alpha/2 *  Q(theta) )
      ul <- -cutoff
      lr <- cutoff
      ur <- theta + qnorm( 1 - alpha/2 *  Q(theta) )
    }
    
    if (theta1 <= theta && theta < theta2) {
      ll <- -cutoff
      ul <- -cutoff
      lr <- cutoff
      ur <- theta + qnorm( pnorm( cutoff - theta ) + (1 - alpha) *  Q(theta) )
    }
    
    if (theta2 <= theta) {
      ll <- -cutoff
      ul <- -cutoff
      lr <- theta - qnorm( .5 * ( 1 + (1 - alpha) *  Q(theta) ) )
      ur <- theta + qnorm( .5 * ( 1 + (1 - alpha) *  Q(theta) ) )
    }
    l <- (ul - ll) + (ur - lr)  #compute the length of the AR
    if(ll==-cutoff) {
      ll <- NA
      ul <- NA
    }
    
    A <- c(ll,ul,lr,ur)
    
  }
  v <- list(A = A, length= l)
  return(v)
}
## Testing:
# ShortestAR(theta=1, cutoff=1, alpha=0.05)




# inputs are: x (value of the observation, |x|>cutoff), 
# sigsq (a known value for the variance of X), 
# cutoff (positive threshold) and 
# alpha (level of the test)
ShortestCI <- function(x, sigsq, cutoff, alpha) {
  stopifnot(abs(x) > cutoff)
  
  if (x<0) {
    tmp <- Recall(x=-x, sigsq=sigsq, cutoff=cutoff, alpha=alpha)
    lower <- -tmp$upper
    upper <- -tmp$lower
  } else {
    
    # reduce the problem to the canonical form: sigsq=1
    x <- x/sqrt(sigsq)
    cutoff <- cutoff/sqrt(sigsq)
    
    Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)	
    #compute useful quantities
    f <- function(theta) ( pnorm(cutoff + theta) - pnorm(cutoff - theta) ) - (1 - alpha) * Q(theta)
    theta1 <- uniroot(f,c(0, cutoff + qnorm(1 - alpha))) $ root
    f <- function(theta) 2*pnorm(theta - cutoff) -1 - (1 - alpha)* Q(theta)
    theta2 <- uniroot( f,c( theta1,cutoff + qnorm(1 - alpha/2) ) ) $ root      # theta2 - cutoff cannot be greater than qnorm(1 - alpha/2)!
    x0 <- qnorm( 1 - .5*alpha * Q(0) )  # not used here, but leaving it in
    x1 <-  theta1 + qnorm ( pnorm(cutoff - theta1) + (1 - alpha)*( 1 - pnorm(cutoff + theta1) + 1 - pnorm(cutoff - theta1) ) )
    x2 <- theta2 + (theta2 - cutoff)
    f <- function(theta) 1 + ( -dnorm(cutoff - theta) + (1 - alpha) * ( -dnorm(cutoff + theta) + dnorm(cutoff - theta) ) ) / ( dnorm( qnorm( pnorm(cutoff - theta) + (1 - alpha) * (1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)) ) ) )
    R <- uniroot(f, c(theta1, theta2)) $ root
    
    #obtain lower end of CI
    if (cutoff < x && x < x1) {
      f <- function(theta) 2 * (1 - pnorm(x - theta)) - alpha *  Q(theta)
      lower <- uniroot(f, c(-theta1 - .1, theta1)) $ root
    }
    if (x1 < x && x < x2) {
      f <- function(theta) pnorm( x - theta ) - pnorm( cutoff - theta ) - (1 - alpha) *  Q(theta)
      lower <- uniroot(f, c(R, theta2)) $ root
    }
    if (x > x2) {
      f <- function(theta) 2 * pnorm( x - theta ) - 1 - (1 - alpha) *  Q(theta)
      lower <- uniroot(f, c(theta2, x)) $ root
    }
    #obtain upper end of CI
    f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) *  Q(theta)
    upper <- uniroot( f, c( theta2, x + 2 * qnorm(1 - alpha/2) ) ) $ root # upper is < x + 1.96, but for root-finding issues I put 2 * 1.96
    
    # rescale
    lower <- sqrt(sigsq)*lower
    upper <- sqrt(sigsq)*upper
    
  }
  CI <- list(lower = lower, upper = upper)
  return(CI)
}
## Testing:
# ShortestCI(x=1.9, sigsq=1, cutoff=1, alpha=0.01) 




# inputs are: theta (value of location parameter), r (the maximum of the ration between the MP interval and the usual two-sided interavl), cutoff (positive threshold) and alpha (level of the test)
ModifiedPrattAR <- function(theta,ratio,cutoff,alpha){
  Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)
  #compute useful quantities
  f <- function(theta) ( pnorm(cutoff + theta) - pnorm(cutoff - theta) ) - (1 - alpha) * Q(theta)
  theta1 <- uniroot(f,c(0,cutoff + qnorm(1 - alpha))) $ root
  f <- function(theta) pnorm( cutoff + ratio * ShortestAR(theta,cutoff,alpha)$l - theta ) - pnorm( cutoff - theta ) - (1 - alpha) *  Q(theta)
  thetatilde1 <- uniroot(f,c(0,theta1)) $ root
  thetatilde2 <- uniroot(f,c(theta1, ratio * 2 * qnorm(1 - alpha / 2)))$root
  
  if (theta<0) {
    tmp <- ModifiedPrattAR(-theta,ratio,cutoff,alpha)
    A <- -rev(tmp)
  } else{
    
    #compute ends of the AR
    if (theta == 0) {
      ll <- -qnorm( 1 - alpha/2 * Q(0) )
      ul <- -cutoff
      lr <- cutoff
      ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
    if (0 < theta && theta < thetatilde1) {
      l <- ShortestAR(theta,cutoff,alpha)$l
      f <- function(x) 1 - pnorm(theta - x) + 1 - pnorm(cutoff + ratio * l - (-cutoff - x) - theta) - alpha *  Q(theta)
      atilde1 <- uniroot(f,c(theta - qnorm(1 - alpha/2 *  Q(theta)), -cutoff))$root
      ll <- atilde1
      ul <- -cutoff
      lr <- cutoff
      ur <- cutoff + ratio * l - (-cutoff - atilde1)
    }
    if (thetatilde1 <= theta && theta < thetatilde2) {
      l <- ShortestAR(theta,cutoff,alpha)$l
      f <- function(x) pnorm(x + ratio * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
      atilde2 <- uniroot(f,c(cutoff, theta - qnorm( (1 - alpha) *  Q(theta) )))$root
      ll <- NA
      ul <- NA
      lr <- atilde2
      ur <- atilde2 + ratio * l
    }
    if (thetatilde2 <= theta) {
      l <- ShortestAR(theta,cutoff,alpha)$l
      f <- function(x) pnorm(x + ratio * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
      atilde2 <- uniroot(f,c(theta - ratio * l /2, theta - qnorm((1 - alpha) *  Q(theta))))$root
      ll <- NA
      ul <- NA
      lr <- atilde2
      ur <- atilde2 + ratio * l
    }
    A <- c(ll,ul,lr,ur)
  }
  return(A)
}
## Testing:
# ModifiedPrattAR(theta=1 ,r ,cutoff, alpha)
  




# inputs are: x (value of observation, |x|>cutoff), sigsq (a known value for the variance of X), r (the maximum of the ration between the MP interval and the usual two-sided interavl), cutoff (positive threshold) and alpha (level of the test)
ModifiedPrattCI <- function(x, sigsq , ratio, cutoff, alpha) {
  
  if (x<0) {
    tmp <- ModifiedPrattCI(-x,sigsq,ratio,cutoff,alpha)
    lower <- -tmp$upper
    upper <- -tmp$lower
  } else {
    
    # reduce the problem to the canonical form: sigsq=1
    x <- x/sqrt(sigsq)
    cutoff <- cutoff/sqrt(sigsq)	
    
    Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)	
    #compute useful quantities
    f <- function(theta) ( pnorm(cutoff + theta) - pnorm(cutoff - theta) ) - (1 - alpha) * Q(theta)
    theta1 <- uniroot(f,c(0,cutoff + qnorm(1 - alpha))) $ root
    f <- function(theta) pnorm( cutoff + ratio * ShortestAR(theta,cutoff,alpha)$l - theta ) - pnorm( cutoff - theta ) - (1 - alpha) *  Q(theta)
    thetatilde1 <- uniroot(f,c(0,theta1)) $ root
    thetatilde2 <- uniroot(f,c(theta1, ratio * 2 * qnorm(1 - alpha / 2)))$root
    zalphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
    lzero <- ShortestAR(0,cutoff,alpha)$l
    f <- function(x) 1 - pnorm(x) + 1 - pnorm(2 * cutoff + ratio * lzero - x) - 2 * alpha * (1 - pnorm(cutoff))
    xtilde1 <- uniroot(f,c(cutoff,zalphahalf))$root
    xtilde2 <- uniroot(f,c(zalphahalf,cutoff + ratio * 2 * qnorm(1 - alpha / 2)))$root
    ltilde1 <- ShortestAR(thetatilde1,cutoff,alpha)$l
    
    #obtain CI ends
    is.neg <- 0
    if (x < 0) is.neg <- 1
    x <- abs(x)
    
    #obtain lower end of CI
    if (cutoff < x && x < xtilde1) {
      f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * cutoff + ratio * ShortestAR(theta,cutoff,alpha)$l) - alpha *  Q(theta)
      lower <- uniroot(f,c(-thetatilde1 - 1e-3,0))$ root      #the -1 in the lower boundary of the interval is for technical reasons (the root really lies in (-thetatilde1,0)  #verify that there exist only one root
    }
    if (xtilde1 <= x && x < zalphahalf) {
      lower <- 0
    }
    if (zalphahalf <= x && x < xtilde2) {
      lower <- 0
    }
    if (xtilde2 <= x && x < cutoff + ratio * ltilde1 ) {
      f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * cutoff + ratio * ShortestAR(theta,cutoff,alpha)$l) - alpha *  Q(theta)
      lower <- uniroot(f,c(0,thetatilde1))$ root
    }
    if (x > cutoff + ratio * ltilde1 ) {
      f <- function(theta) pnorm(x - theta) - pnorm(x - ratio * ShortestAR(theta,cutoff,alpha)$l - theta) - (1 - alpha) *  Q(theta)
      m <- optimize(f, c(thetatilde1,x), maximum=T)$maximum
      lower <- uniroot(f, c(thetatilde1,m))$root
    }
    #obtain upper end of CI
    f <- function(theta) pnorm(x + ratio * ShortestAR(theta,cutoff,alpha)$l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
    m <- optimize(f, c(x, x + ratio * 2 * qnorm(1 - alpha / 2)), maximum = T)$maximum
    upper <- uniroot(f,c(0,m))$root
    
    # rescale
    lower <- sqrt(sigsq)*lower
    upper <- sqrt(sigsq)*upper
    
  }
  CI <- list(lower = lower, upper = upper)
  return(CI)
}
## Testing:
#ModifiedPrattCI(10, 1 , ratio=1.5, cutoff=1, alpha=0.05)
  




# inputs are: theta (value of location parameter), 
# lambda (positive penalty term for length of the AR), 
# cutoff (positive threshold) and 
# alpha (level of the test)
QuasiConventionalAR <- function(theta, lambda, cutoff, alpha){
  Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)
  #compute useful quantities
  f <- function(theta) ( pnorm(cutoff + theta) - pnorm(cutoff - theta) ) - (1 - alpha) * Q(theta)
  theta1 <- uniroot(f,c(0, cutoff + qnorm(1 - alpha))) $ root
  f <- function(theta) 1 - pnorm(cutoff - theta) - (1 - alpha) * Q(theta)
  thetastar <- uniroot(f,c(0,cutoff + qnorm(1 - alpha)))$root
  f <- function(theta) 1 + lambda * ( 1 - dnorm(cutoff + theta) / dnorm(qnorm(2 - pnorm(cutoff + theta) - alpha * Q(theta))) )
  b <- thetastar
  while(is.na(f(b))==T) {b <- b + 1e-5}
  thetaprime1 <- uniroot(f, c(b, theta1))$root
  
  if (theta<0) {
    tmp <- Recall(theta=-theta, lambda=lambda, cutoff=cutoff, alpha=alpha)
    A <- -rev(tmp)
  } else {
    
    #compute ends of the AR
    
    if (theta == 0) {
      ll <- -qnorm( 1 - alpha/2 * Q(0) )
      ul <- -cutoff
      lr <- cutoff
      ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
    if (0 < theta && theta < thetaprime1) {
      f <- function(d) 1 + lambda * ( 1 - dnorm(cutoff + d + theta) / dnorm(qnorm(2 - pnorm(cutoff + d + theta) - alpha * Q(theta))) )
      dunderbar <- max(-cutoff - theta + qnorm(1 - alpha * Q(theta)), 0)
      dbar <- -cutoff - ShortestAR(theta,cutoff,alpha)$A[1]
      dstar <- uniroot(f,c(dunderbar+1e-5, dbar))$root
      aprime1 <- theta - cutoff + qnorm( 1 - pnorm(cutoff + dstar + theta) + 1 - alpha * Q(theta) )
      ll <- -cutoff - dstar
      ul <- -cutoff
      lr <- cutoff
      ur <- cutoff + aprime1
    }
    if (thetaprime1 <= theta && theta < theta1) {
      aprime2 <- theta - cutoff + qnorm( pnorm(cutoff - theta) + (1 - alpha) * Q(theta) )
      ll <- NA
      ul <- NA
      lr <- cutoff
      ur <- cutoff + aprime2
    }
    if (theta > theta1) {
      A <- ShortestAR(theta,cutoff,alpha)$A
      ll <- A[1]
      ul <- A[2]
      lr <- A[3]
      ur <- A[4]
    }
    A <- c(ll, ul, lr, ur)
    
  }
  return(A)
}




# inputs are: x (value of observation, |x|>cutoff), 
# sigsq (a known value for the variance of X), 
# lambada (positive pentalty term for length of the AR), 
# cutoff (positive threshold) and 
# alpha (level of the test)
QuasiConventionalCI <- function(x, sigsq, lambda ,cutoff ,alpha){
  if (x<0) {
    tmp <- Recall(x=-x, sigsq=sigsq, lambda=lambda, cutoff=cutoff, alpha=alpha)
    lower <- -tmp$upper
    upper <- -tmp$lower
  } else {
    
    # reduce the problem to the canonical form: sigsq=1
    x <- x/sqrt(sigsq)
    cutoff <- cutoff/sqrt(sigsq)	
    
    Q <- function(theta) 1 - pnorm(cutoff + theta) + 1 - pnorm(cutoff - theta)
    f <- function(theta) ( pnorm(cutoff + theta) - pnorm(cutoff - theta) ) - (1 - alpha) * Q(theta)
    theta1 <- uniroot(f,c(0, cutoff + qnorm(1 - alpha))) $ root
    f <- function(theta) 1 - pnorm(cutoff - theta) - (1 - alpha) * Q(theta)
    thetastar <- uniroot(f,c(0,cutoff + qnorm(1 - alpha)))$root
    f <- function(theta) 1 + lambda * ( 1 - dnorm(cutoff + theta) / dnorm(qnorm(2 - pnorm(cutoff + theta) - alpha * Q(theta))) )
    b <- thetastar
    while(is.na(f(b))==T) {b <- b + 1e-5} # here it would be nice to suppress warnings(AW)
    thetaprime1 <- uniroot(f, c(b, theta1))$root
    
    dmin <- function(theta) max( -cutoff - theta + qnorm(1 - alpha * Q(theta)), 0 )
    dmax <- function(theta) -cutoff - ShortestAR(theta,cutoff,alpha)$A[1]
    #compute useful quantities
    zalphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
    dminzero <- dmin(0)
    dmaxzero <- dmax(0)
    f <- function(d) 1 + lambda * ( 1 - dnorm(cutoff + d) / dnorm(qnorm(2 - pnorm(cutoff + d) - alpha * Q(0))) )
    xprime1 <- uniroot(f, c(dminzero, dmaxzero))$root + cutoff
    xprime2 <- qnorm( 2 - pnorm(xprime1) - alpha * Q(0) )
    xprime3 <- thetaprime1 + qnorm( pnorm(cutoff - thetaprime1) + (1 - alpha) * Q(thetaprime1) )
    
    #obtain CI ends
    is.neg <- 0
    if (x < 0) is.neg <- 1
    x <- abs(x)
    
    #obtain lower end of CI        #NOTE: numerical problem in obtaining lower end when |x-cutoff| extremely small
    if (cutoff < x && x < xprime1) {
      d <- x - cutoff
      f <- function(theta) 1 + lambda * ( 1 - dnorm(cutoff + d + theta) / dnorm(qnorm(2 - pnorm(cutoff + d + theta) - alpha * Q(theta))) )
      g <- function(theta) dmin(theta) - (x - cutoff)
      if (x - cutoff > dmin(0)) {thetamin <- 0} else {
        g <- function(theta) dmin(theta) - (x - cutoff)
        thetamin <- uniroot(g, c(0, thetastar + 1e-3))$root
      }
      leftend <- thetamin + 1e-5
      while (is.na(suppressWarnings(f(leftend)))) leftend <- leftend + 1e-5 * 5
      lower <- uniroot(f,c(leftend, thetaprime1 + .1)) $ root       #NOTE: for very small x, this might result in wrong value of theta for the lower bound - one that is smaller than -thetaprime1 !
      #(but I am willing to live with that.. this value is in any case very close to the correct one,  which is about -thetaprime1)
      lower <- - lower
    }
    if (xprime1 <= x && x < zalphahalf) {
      lower <- 0    # 0 is included in CI
    }
    if (zalphahalf <= x && x < xprime2) {
      lower <- 0    # 0 is NOT included in CI
    }
    if (xprime2 <= x && x < xprime3) {
      g <- function(theta) -cutoff - theta + qnorm( 2 - pnorm(x - theta) - alpha * Q(theta) )
      f <- function(theta) 1 + lambda * ( 1 - dnorm(cutoff + g(theta) + theta) / dnorm(qnorm(2 - pnorm(cutoff + g(theta) + theta) - alpha * Q(theta))) )
      h <- function(theta) 1 - pnorm(x - theta) - alpha * Q(theta)
      thetamax <- uniroot(h,c(0,x))$root
      lower <- uniroot(f,c(0 - 1e-1, thetamax - 1e-1)) $ root
    }
    if (x >= xprime3 ) {
      lower <- ShortestCI(x=x, sigsq=sigsq, cutoff=cutoff, alpha=alpha)$lower
    }
    
    #obtain upper end of CI
    f <- function(theta) 2*pnorm(theta - cutoff) -1 - (1 - alpha)* Q(theta)
    theta2 <- uniroot( f,c( theta1,cutoff + qnorm(1 - alpha/2) ) ) $ root      # theta2 - cutoff cannot be greater than qnorm(1 - alpha/2)!
    f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) * Q(theta)
    upper <- uniroot( f, c( theta2, x + 1.1 * qnorm(1 - alpha/2) ) ) $ root # upper is < x + 1.96, but for root-finding issues I put 2 * 1.96
    
    # rescale
    lower <- sqrt(sigsq)*lower
    upper <- sqrt(sigsq)*upper
    
  }
  CI <- list(lower = lower, upper = upper)
  return(CI)
}
## Testing:
#QuasiConventionalCI(x=-2, sigsq=1, lambda=2, cutoff=1, alpha=0.05)









# the following functions implement methods suggested by Zhong and Prentice (2008). 
# Each directly computes bounds for a conditional CI (without ARs)
# Computed only for the case sigsq=1

# this one is based on the likelihood ratio
LikelihoodRatioCI <- function(x, cutoff, alpha){
  log.cond.lik <- function(theta) log( dnorm(x - theta) / (1 - pnorm(cutoff - theta) + 1 - pnorm(cutoff + theta)) )
  g <- function(theta) (2 * pi) ^ (-1/2) * (x - theta) * exp(-(x - theta) ^ 2 / 2) * (1 - pnorm(cutoff - theta) + 
                                                                                        1 - pnorm(cutoff + theta)) - dnorm(x - theta) * (dnorm(cutoff - theta) - dnorm(cutoff + theta))
  root <- uniroot(g,c(min(x,0), max(x,0))) $ root
  h <- function(theta) log.cond.lik(theta) - ( log.cond.lik(root) - qgamma(1 - alpha, shape=1/2, scale=2) / 2 )
  
  #this is a technical part i'm using because I can't specify lower end for root finding as +/-inf
  t1 <- root
  while (h(t1) > 0) t1 <- t1 - 1
  t2 <- root
  while (h(t2) > 0) t2 <- t2 + 1
  
  lower <- uniroot(h, c(t1,root)) $ root
  upper <- uniroot(h, c(root,t2)) $ root
  CI <- list(lower = lower, upper = upper)
  return(CI)
}




# this CI is quantile-based
QuantileCI <- function(x, cutoff, alpha){
  f <- function(theta) 1 - pnorm(x - theta) - alpha / 2 * Q(theta)
  lower <- uniroot(f, c(-cutoff - qnorm(1 - alpha / 2), x - qnorm(1 - alpha / 2))) $ root
  f <- function(theta) pnorm(theta - x) - (1 - alpha / 2) * Q(theta)
  upper <- uniroot(f, c(x, x + qnorm(1 - alpha / 2))) $ root
  CI <- list(lower = lower, upper = upper)
  return(CI)
}





