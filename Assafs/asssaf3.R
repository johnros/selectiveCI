#CCI functions
#commentary added Aug 16th, 2013.

  # this functions computes P(|X|>c) as a function of theta. It is called by functions below
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)

  # The following functions implement the three methods suggested in the paper to construct conditional CIs for theta. The first method is "Shortest-AR", the second is "Modified Pratt" and the third is "Quasi-Conventional". There are two functions corresponding to each of the three methods: The first computes the acceptance regions (AR), and the second computes the CI (that are obtained when inverting the ARs and taking the convex hull).


  # inputs are: theta (value of location parameter), c (positive threshold) and alpha (level of the test)
  # Note: ARs are constructed for sigsq=1; This is enough to obtain a CI for a general sigsq (proper modification is made in the CI function)
Shortest.AR <- function(theta,c,alpha){
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)	
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0,c + qnorm(1 - alpha))) $ root
f <- function(theta) 2*pnorm(theta - c) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,c + qnorm(1 - alpha/2) ) ) $ root      # theta2 - c cannot be greater than qnorm(1 - alpha/2)!

if (theta<0) {
	tmp <- Shortest.AR(-theta,c,alpha)
	A <- -rev(tmp$A)
	l <- tmp$l
} else{

#compute ends of the AR
  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -c
  lr <- c
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }

  if (0 < theta && theta < theta1) {
  ll <- theta - qnorm( 1 - alpha/2 *  Q(theta) )
  ul <- -c
  lr <- c
  ur <- theta + qnorm( 1 - alpha/2 *  Q(theta) )
    }

  if (theta1 <= theta && theta < theta2) {
  ll <- -c
  ul <- -c
  lr <- c
  ur <- theta + qnorm( pnorm( c - theta ) + (1 - alpha) *  Q(theta) )
    }

  if (theta2 <= theta) {
  ll <- -c
  ul <- -c
  lr <- theta - qnorm( .5 * ( 1 + (1 - alpha) *  Q(theta) ) )
  ur <- theta + qnorm( .5 * ( 1 + (1 - alpha) *  Q(theta) ) )
    }
l <- (ul - ll) + (ur - lr)  #compute the length of the AR
if(ll==-c) {
	ll <- NA
	ul <- NA
}

A <- c(ll,ul,lr,ur)

}
v <- list(A = A,l = l)
return(v)
}




  # inputs are: x (value of the observation, |x|>c), sigsq (a known value for the variance of X), c (positive threshold) and alpha (level of the test)
Shortest.CI <- function(x,sigsq,c,alpha) {

if (x<0) {
	tmp <- Shortest.CI(-x,sigsq,c,alpha)
	lower <- -tmp$upper
	upper <- -tmp$lower
} else {

# reduce the problem to the canonical form: sigsq=1
x <- x/sqrt(sigsq)
c <- c/sqrt(sigsq)

Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)	
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, c + qnorm(1 - alpha))) $ root
f <- function(theta) 2*pnorm(theta - c) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,c + qnorm(1 - alpha/2) ) ) $ root      # theta2 - c cannot be greater than qnorm(1 - alpha/2)!
x0 <- qnorm( 1 - .5*alpha * Q(0) )  # not used here, but leaving it in
x1 <-  theta1 + qnorm ( pnorm(c - theta1) + (1 - alpha)*( 1 - pnorm(c + theta1) + 1 - pnorm(c - theta1) ) )
x2 <- theta2 + (theta2 - c)
f <- function(theta) 1 + ( -dnorm(c - theta) + (1 - alpha) * ( -dnorm(c + theta) + dnorm(c - theta) ) ) / ( dnorm( qnorm( pnorm(c - theta) + (1 - alpha) * (1 - pnorm(c + theta) + 1 - pnorm(c - theta)) ) ) )
R <- uniroot(f, c(theta1, theta2)) $ root
	
#obtain lower end of CI
if (c < x && x < x1) {
  f <- function(theta) 2 * (1 - pnorm(x - theta)) - alpha *  Q(theta)
  lower <- uniroot(f, c(-theta1 - .1, theta1)) $ root
  }
if (x1 < x && x < x2) {
  f <- function(theta) pnorm( x - theta ) - pnorm( c - theta ) - (1 - alpha) *  Q(theta)
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




  # inputs are: theta (value of location parameter), r (the maximum of the ration between the MP interval and the usual two-sided interavl), c (positive threshold) and alpha (level of the test)
MP.AR <- function(theta,r,c,alpha){
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0,c + qnorm(1 - alpha))) $ root
f <- function(theta) pnorm( c + r * Shortest.AR(theta,c,alpha)$l - theta ) - pnorm( c - theta ) - (1 - alpha) *  Q(theta)
thetatilde1 <- uniroot(f,c(0,theta1)) $ root
thetatilde2 <- uniroot(f,c(theta1, r * 2 * qnorm(1 - alpha / 2)))$root

if (theta<0) {
	tmp <- MP.AR(-theta,r,c,alpha)
	A <- -rev(tmp)
} else{

#compute ends of the AR
  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -c
  lr <- c
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
  if (0 < theta && theta < thetatilde1) {
l <- Shortest.AR(theta,c,alpha)$l
f <- function(x) 1 - pnorm(theta - x) + 1 - pnorm(c + r * l - (-c - x) - theta) - alpha *  Q(theta)
atilde1 <- uniroot(f,c(theta - qnorm(1 - alpha/2 *  Q(theta)), -c))$root
  ll <- atilde1
  ul <- -c
  lr <- c
  ur <- c + r * l - (-c - atilde1)
    }
  if (thetatilde1 <= theta && theta < thetatilde2) {
l <- Shortest.AR(theta,c,alpha)$l
f <- function(x) pnorm(x + r * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
atilde2 <- uniroot(f,c(c, theta - qnorm( (1 - alpha) *  Q(theta) )))$root
  ll <- NA
  ul <- NA
  lr <- atilde2
  ur <- atilde2 + r * l
    }
  if (thetatilde2 <= theta) {
l <- Shortest.AR(theta,c,alpha)$l
f <- function(x) pnorm(x + r * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
atilde2 <- uniroot(f,c(theta - r * l /2, theta - qnorm((1 - alpha) *  Q(theta))))$root
  ll <- NA
  ul <- NA
  lr <- atilde2
  ur <- atilde2 + r * l
    }
A <- c(ll,ul,lr,ur)
}
return(A)
}




  # inputs are: x (value of observation, |x|>c), sigsq (a known value for the variance of X), r (the maximum of the ration between the MP interval and the usual two-sided interavl), c (positive threshold) and alpha (level of the test)
MP.CI <- function(x,sigsq,r,c,alpha) {

if (x<0) {
	tmp <- MP.CI(-x,sigsq,r,c,alpha)
	lower <- -tmp$upper
	upper <- -tmp$lower
} else {

# reduce the problem to the canonical form: sigsq=1
x <- x/sqrt(sigsq)
c <- c/sqrt(sigsq)	
	
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)	
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0,c + qnorm(1 - alpha))) $ root
f <- function(theta) pnorm( c + r * Shortest.AR(theta,c,alpha)$l - theta ) - pnorm( c - theta ) - (1 - alpha) *  Q(theta)
thetatilde1 <- uniroot(f,c(0,theta1)) $ root
thetatilde2 <- uniroot(f,c(theta1, r * 2 * qnorm(1 - alpha / 2)))$root
zalphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
lzero <- Shortest.AR(0,c,alpha)$l
f <- function(x) 1 - pnorm(x) + 1 - pnorm(2 * c + r * lzero - x) - 2 * alpha * (1 - pnorm(c))
xtilde1 <- uniroot(f,c(c,zalphahalf))$root
xtilde2 <- uniroot(f,c(zalphahalf,c + r * 2 * qnorm(1 - alpha / 2)))$root
ltilde1 <- Shortest.AR(thetatilde1,c,alpha)$l

#obtain CI ends
is.neg <- 0
  if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI
if (c < x && x < xtilde1) {
  f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * c + r * Shortest.AR(theta,c,alpha)$l) - alpha *  Q(theta)
  lower <- uniroot(f,c(-thetatilde1 - 1e-3,0))$ root      #the -1 in the lower boundary of the interval is for technical reasons (the root really lies in (-thetatilde1,0)  #verify that there exist only one root
  }
if (xtilde1 <= x && x < zalphahalf) {
  lower <- 0
  }
if (zalphahalf <= x && x < xtilde2) {
  lower <- 0
  }
if (xtilde2 <= x && x < c + r * ltilde1 ) {
  f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * c + r * Shortest.AR(theta,c,alpha)$l) - alpha *  Q(theta)
  lower <- uniroot(f,c(0,thetatilde1))$ root
  }
if (x > c + r * ltilde1 ) {
  f <- function(theta) pnorm(x - theta) - pnorm(x - r * Shortest.AR(theta,c,alpha)$l - theta) - (1 - alpha) *  Q(theta)
  m <- optimize(f, c(thetatilde1,x), maximum=T)$maximum
  lower <- uniroot(f, c(thetatilde1,m))$root
  }
#obtain upper end of CI
f <- function(theta) pnorm(x + r * Shortest.AR(theta,c,alpha)$l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
m <- optimize(f, c(x, x + r * 2 * qnorm(1 - alpha / 2)), maximum = T)$maximum
upper <- uniroot(f,c(0,m))$root

# rescale
lower <- sqrt(sigsq)*lower
upper <- sqrt(sigsq)*upper

}
CI <- list(lower = lower, upper = upper)
return(CI)
}




  # inputs are: theta (value of location parameter), lambad (positive penalty term for length of the AR), c (positive threshold) and alpha (level of the test)
QC.AR <- function(theta,lambda,c,alpha){
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, c + qnorm(1 - alpha))) $ root
f <- function(theta) 1 - pnorm(c - theta) - (1 - alpha) * Q(theta)
thetastar <- uniroot(f,c(0,c + qnorm(1 - alpha)))$root
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + theta) / dnorm(qnorm(2 - pnorm(c + theta) - alpha * Q(theta))) )
b <- thetastar
while(is.na(f(b))==T) {b <- b + 1e-5}
thetaprime1 <- uniroot(f, c(b, theta1))$root

if (theta<0) {
	tmp <- QC.AR(-theta,lambda,c,alpha)
	A <- -rev(tmp)
} else {

#compute ends of the AR

  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -c
  lr <- c
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
  if (0 < theta && theta < thetaprime1) {
f <- function(d) 1 + lambda * ( 1 - dnorm(c + d + theta) / dnorm(qnorm(2 - pnorm(c + d + theta) - alpha * Q(theta))) )
dunderbar <- max(-c - theta + qnorm(1 - alpha * Q(theta)), 0)
dbar <- -c - Shortest.AR(theta,c,alpha)$A[1]
dstar <- uniroot(f,c(dunderbar+1e-5, dbar))$root
aprime1 <- theta - c + qnorm( 1 - pnorm(c + dstar + theta) + 1 - alpha * Q(theta) )
  ll <- -c - dstar
  ul <- -c
  lr <- c
  ur <- c + aprime1
    }
  if (thetaprime1 <= theta && theta < theta1) {
aprime2 <- theta - c + qnorm( pnorm(c - theta) + (1 - alpha) * Q(theta) )
  ll <- NA
  ul <- NA
  lr <- c
  ur <- c + aprime2
    }
  if (theta > theta1) {
A <- Shortest.AR(theta,c,alpha)$A
  ll <- A[1]
  ul <- A[2]
  lr <- A[3]
  ur <- A[4]
    }
A <- c(ll, ul, lr, ur)

}
return(A)
}




  # inputs are: x (value of observation, |x|>c), sigsq (a known value for the variance of X), lambada (positive pentalty term for length of the AR), c (positive threshold) and alpha (level of the test)
QC.CI <- function(x,sigsq,lambda,c,alpha){

if (x<0) {
	tmp <- QC.CI(-x,sigsq,lambda,c,alpha)
	lower <- -tmp$upper
	upper <- -tmp$lower
} else {

# reduce the problem to the canonical form: sigsq=1
x <- x/sqrt(sigsq)
c <- c/sqrt(sigsq)	
	
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, c + qnorm(1 - alpha))) $ root
f <- function(theta) 1 - pnorm(c - theta) - (1 - alpha) * Q(theta)
thetastar <- uniroot(f,c(0,c + qnorm(1 - alpha)))$root
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + theta) / dnorm(qnorm(2 - pnorm(c + theta) - alpha * Q(theta))) )
b <- thetastar
while(is.na(f(b))==T) {b <- b + 1e-5} # here it would be nice to suppress warnings(AW)
thetaprime1 <- uniroot(f, c(b, theta1))$root

dmin <- function(theta) max( -c - theta + qnorm(1 - alpha * Q(theta)), 0 )
dmax <- function(theta) -c - Shortest.AR(theta,c,alpha)$A[1]
#compute useful quantities
zalphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
dminzero <- dmin(0)
dmaxzero <- dmax(0)
f <- function(d) 1 + lambda * ( 1 - dnorm(c + d) / dnorm(qnorm(2 - pnorm(c + d) - alpha * Q(0))) )
xprime1 <- uniroot(f, c(dminzero, dmaxzero))$root + c
xprime2 <- qnorm( 2 - pnorm(xprime1) - alpha * Q(0) )
xprime3 <- thetaprime1 + qnorm( pnorm(c - thetaprime1) + (1 - alpha) * Q(thetaprime1) )

#obtain CI ends
is.neg <- 0
if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI        #NOTE: numerical problem in obtaining lower end when |x-c| extremely small
if (c < x && x < xprime1) {
d <- x - c
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + d + theta) / dnorm(qnorm(2 - pnorm(c + d + theta) - alpha * Q(theta))) )
g <- function(theta) dmin(theta) - (x - c)
  if (x - c > dmin(0)) {thetamin <- 0} else {
    g <- function(theta) dmin(theta) - (x - c)
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
g <- function(theta) -c - theta + qnorm( 2 - pnorm(x - theta) - alpha * Q(theta) )
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + g(theta) + theta) / dnorm(qnorm(2 - pnorm(c + g(theta) + theta) - alpha * Q(theta))) )
h <- function(theta) 1 - pnorm(x - theta) - alpha * Q(theta)
thetamax <- uniroot(h,c(0,x))$root
lower <- uniroot(f,c(0 - 1e-1, thetamax - 1e-1)) $ root
  }
if (x >= xprime3 ) {
  lower <- Shortest.CI(x,c,alpha)$lower
  }
  
#obtain upper end of CI
f <- function(theta) 2*pnorm(theta - c) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,c + qnorm(1 - alpha/2) ) ) $ root      # theta2 - c cannot be greater than qnorm(1 - alpha/2)!
f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) * Q(theta)
upper <- uniroot( f, c( theta2, x + 1.1 * qnorm(1 - alpha/2) ) ) $ root # upper is < x + 1.96, but for root-finding issues I put 2 * 1.96

# rescale
lower <- sqrt(sigsq)*lower
upper <- sqrt(sigsq)*upper

}
CI <- list(lower = lower, upper = upper)
return(CI)
}









  # the following functions implement methods suggested by Zhong and Prentice (2008). Each directly computes bounds for a conditional CI (without ARs)
  # Computed only for the case sigsq=1

  # this one is based on the likelihood ratio
LR.CI <- function(x,c,alpha){
log.cond.lik <- function(theta) log( dnorm(x - theta) / (1 - pnorm(c - theta) + 1 - pnorm(c + theta)) )
g <- function(theta) (2 * pi) ^ (-1/2) * (x - theta) * exp(-(x - theta) ^ 2 / 2) * (1 - pnorm(c - theta) + 
1 - pnorm(c + theta)) - dnorm(x - theta) * (dnorm(c - theta) - dnorm(c + theta))
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
QB.CI <- function(x,c,alpha){

f <- function(theta) 1 - pnorm(x - theta) - alpha / 2 * Q(theta)
lower <- uniroot(f, c(-c - qnorm(1 - alpha / 2), x - qnorm(1 - alpha / 2))) $ root
f <- function(theta) pnorm(theta - x) - (1 - alpha / 2) * Q(theta)
upper <- uniroot(f, c(x, x + qnorm(1 - alpha / 2))) $ root
CI <- list(lower = lower, upper = upper)
return(CI)
}



  # ignore this part

# Length of the QC CI at x = xprime2 VS lambda		#<-- garbage?

QC.maxlength <- function(lambda,c,alpha){
Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)
dmin <- function(theta) max( -c - theta + qnorm(1 - alpha * Q(theta)), 0 )
dmax <- function(theta) -c - Shortest.AR(theta,c,alpha)$A[1]
f <- function(d) 1 + lambda * ( 1 - dnorm(c + d) / dnorm(qnorm(2 - pnorm(c + d) - alpha * Q(0))) )
dminzero <- dmin(0)
dmaxzero <- dmax(0)
while (is.na(suppressWarnings(f(dminzero)))) dminzero <- dminzero + 1e-12  #just a technical issue with computing f() exactly at dminzero	
xprime1 <- uniroot(f, c(dminzero, dmaxzero))$root + c

xprime2 <- qnorm( 2 - pnorm(xprime1) - alpha * Q(0) )
a <- QC.CI(xprime2,lambda,c,alpha) $ lower
b <- QC.CI(xprime2,lambda,c,alpha) $ upper
l <- b - a
return(l)
}



# #   # Examples:
# alpha <- .05
# c <- 1

# # QC method:

# lambda <- 2
# QC.AR(.5,lambda=2,c,alpha) # The AR is a disjoint union of an interval lying to the left of -c and an interval lying to the right of c. output lists lower left, upper left, lower right, upper right endpoints, respectively.
# QC.CI(1.7, sigsq=4,lambda=2,c,alpha)


