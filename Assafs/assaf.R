#functions

Q <- function(theta) 1 - pnorm(c + theta) + 1 - pnorm(c - theta)

Shortest.AR <- function(theta,c,alpha){
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0,c + qnorm(1 - alpha))) $ root
f <- function(theta) 2*pnorm(theta - c) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,c + qnorm(1 - alpha/2) ) ) $ root      # theta2 - c cannot be greater than qnorm(1 - alpha/2)!

#compute ends of the AR
is.neg <- 0
  if (theta < 0) is.neg <- 1
theta <- abs(theta)

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
A <- c(ll,ul,lr,ur)
if (is.neg == 1) A <- c(-ur,-lr,-ul,-ll)
l <- (ul - ll) + (ur - lr)  #compute the length of the AR
v <- list(A = A,l = l)
return(v)
}





Shortest.CI <- function(x,c,alpha) {
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

#obtain CI ends
is.neg <- 0
  if (x < 0) is.neg <- 1
x <- abs(x)

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
CI <- list(lower = lower, upper = upper)
if (is.neg == 1) CI <- list(lower = -upper, upper = -lower)
return(CI)
}





MP.AR <- function(theta,r,c,alpha){
#compute useful quantities
f <- function(theta) pnorm( c + r * Shortest.AR(theta,c,alpha)$l - theta ) - pnorm( c - theta ) - (1 - alpha) *  Q(theta)
thetatilde1 <- uniroot(f,c(0,theta1)) $ root
thetatilde2 <- uniroot(f,c(theta1, r * 2 * qnorm(1 - alpha / 2)))$root

#compute ends of the AR
is.neg <- 0
  if (theta < 0) is.neg <- 1
theta <- abs(theta)

  if (theta == 0) {
  ll <- -qnorm( 1 - alpha/2 * Q(0) )
  ul <- -c
  lr <- c
  ur <- qnorm( 1 - alpha/2 * Q(0) )
    }
  if (0 < theta && theta < thetatilde1) {
l <- Shortest.AR(theta,c,alpha)$l
f <- function(x) 1 - pnorm(theta - x) + 1 - pnorm(c + r * l - (-c - x) - theta) - alpha *  Q(theta)
g1 <- uniroot(f,c(theta - qnorm(1 - alpha/2 *  Q(theta)), -c))$root
  ll <- g1
  ul <- -c
  lr <- c
  ur <- c + r * l - (-c - g1)
    }
  if (thetatilde1 <= theta && theta < thetatilde2) {
l <- Shortest.AR(theta,c,alpha)$l
f <- function(x) pnorm(x + r * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
g2 <- uniroot(f,c(c, theta - qnorm( (1 - alpha) *  Q(theta) )))$root
  ll <- -c
  ul <- -c
  lr <- g2
  ur <- g2 + r * l
    }
  if (thetatilde2 <= theta) {
l <- Shortest.AR(theta,c,alpha)$l
f <- function(x) pnorm(x + r * l - theta) - pnorm(x - theta) - (1 - alpha) *  Q(theta)
g2 <- uniroot(f,c(theta - r * l /2, theta - qnorm((1 - alpha) *  Q(theta))))$root
  ll <- -c
  ul <- -c
  lr <- g2
  ur <- g2 + r * l
    }
A <- c(ll,ul,lr,ur)
if (is.neg == 1) A <- c(-ur,-lr,-ul,-ll)
return(A)
}





MP.CI <- function(x,r,c,alpha){
#compute useful quantities
f <- function(theta) pnorm( c + r * Shortest.AR(theta,c,alpha)$l - theta ) - pnorm( c - theta ) - (1 - alpha) *  Q(theta)
thetatilde1 <- uniroot(f,c(0,theta1)) $ root
thetatilde2 <- uniroot(f,c(theta1, r * 2 * qnorm(1 - alpha / 2)))$root
calphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
lzero <- Shortest.AR(0,c,alpha)$l
f <- function(x) 1 - pnorm(x) + 1 - pnorm(2 * c + r * lzero - x) - 2 * alpha * (1 - pnorm(c))
cbarzero <- uniroot(f,c(c,calphahalf))$root
ctildezero <- uniroot(f,c(calphahalf,c + r * 2 * qnorm(1 - alpha / 2)))$root
ltilde1 <- Shortest.AR(thetatilde1,c,alpha)$l

#obtain CI ends
is.neg <- 0
  if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI
if (c < x && x < cbarzero) {
  f <- function(theta) 1 - pnorm(x - theta) + 1 - pnorm(theta - x + 2 * c + r * Shortest.AR(theta,c,alpha)$l) - alpha *  Q(theta)
  lower <- uniroot(f,c(-thetatilde1 - 1e-3,0))$ root      #the -1 in the lower boundary of the interval is for technical reasons (the root really lies in (-thetatilde1,0)  #verify that there exist only one root
  }
if (cbarzero <= x && x < calphahalf) {
  lower <- 0
  }
if (calphahalf <= x && x < ctildezero) {
  lower <- 0
  }
if (ctildezero <= x && x < c + r * ltilde1 ) {
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
CI <- list(lower = lower, upper = upper)
if (is.neg == 1) CI <- list(lower = -upper, upper = -lower)
return(CI)
}





QC.AR <- function(theta,lambda,c,alpha){
#compute useful quantities
f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, c + qnorm(1 - alpha))) $ root
f <- function(theta) 1 - pnorm(c - theta) - (1 - alpha) * Q(theta)
thetastar <- uniroot(f,c(0,c + qnorm(1 - alpha)))$root
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + theta) / dnorm(qnorm(2 - pnorm(c + theta) - alpha * Q(theta))) )
b <- thetastar
while(is.na(f(b))==T) {b <- b + 1e-5}
thetaprime1 <- uniroot(f, c(b, theta1))$root

#compute ends of the AR
is.neg <- 0
  if (theta < 0) is.neg <- 1
theta <- abs(theta)

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
h1 <- theta - c + qnorm( 1 - pnorm(c + dstar + theta) + 1 - alpha * Q(theta) )
  ll <- -c - dstar
  ul <- -c
  lr <- c
  ur <- c + h1
    }
  if (thetaprime1 <= theta && theta < theta1) {
h2 <- theta - c + qnorm( pnorm(c - theta) + (1 - alpha) * Q(theta) )
  ll <- -c
  ul <- -c
  lr <- c
  ur <- c + h2
    }
  if (theta > theta1) {
A <- Shortest.AR(theta,c,alpha)$A
  ll <- A[1]
  ul <- A[2]
  lr <- A[3]
  ur <- A[4]
    }
A <- c(ll, ul, lr, ur)
if (is.neg == 1) A <- c(-ur,-lr,-ul,-ll)
return(A)
}





QC.CI <- function(x,lambda,c,alpha){

f <- function(theta) ( pnorm(c + theta) - pnorm(c - theta) ) - (1 - alpha) * Q(theta)
theta1 <- uniroot(f,c(0, c + qnorm(1 - alpha))) $ root
f <- function(theta) 1 - pnorm(c - theta) - (1 - alpha) * Q(theta)
thetastar <- uniroot(f,c(0,c + qnorm(1 - alpha)))$root
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + theta) / dnorm(qnorm(2 - pnorm(c + theta) - alpha * Q(theta))) )
b <- thetastar
while(is.na(f(b))==T) {b <- b + 1e-5}
thetaprime1 <- uniroot(f, c(b, theta1))$root

dmin <- function(theta) max( -c - theta + qnorm(1 - alpha * Q(theta)), 0 )
dmax <- function(theta) -c - Shortest.AR(theta,c,alpha)$A[1]
#compute useful quantities
calphahalf <- qnorm( 1 - .5 * alpha * Q(0) )
dminzero <- dmin(0)
dmaxzero <- dmax(0)
f <- function(d) 1 + lambda * ( 1 - dnorm(c + d) / dnorm(qnorm(2 - pnorm(c + d) - alpha * Q(0))) )
cbarzero.qc <- uniroot(f, c(dminzero, dmaxzero))$root + c
ctildezero.qc <- qnorm( 2 - pnorm(cbarzero.qc) - alpha * Q(0) )
ctildethetaprime1 <- thetaprime1 + qnorm( pnorm(c - thetaprime1) + (1 - alpha) * Q(thetaprime1) )

#obtain CI ends
is.neg <- 0
if (x < 0) is.neg <- 1
x <- abs(x)

#obtain lower end of CI        #NOTE: numerical problem in obtaining lower end when |x-c| extremely small
if (c < x && x < cbarzero.qc) {
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
if (cbarzero.qc <= x && x < calphahalf) {
  lower <- 0    # 0 is included in CI
  }
if (calphahalf <= x && x < ctildezero.qc) {
  lower <- 0    # 0 is NOT included in CI
  }
if (ctildezero.qc <= x && x < ctildethetaprime1) {
g <- function(theta) -c - theta + qnorm( 2 - pnorm(x - theta) - alpha * Q(theta) )
f <- function(theta) 1 + lambda * ( 1 - dnorm(c + g(theta) + theta) / dnorm(qnorm(2 - pnorm(c + g(theta) + theta) - alpha * Q(theta))) )
h <- function(theta) 1 - pnorm(x - theta) - alpha * Q(theta)
thetamax <- uniroot(h,c(0,x))$root
lower <- uniroot(f,c(0 - 1e-2, thetamax - 1e-2)) $ root
  }
if (x >= ctildethetaprime1 ) {
  lower <- Shortest.CI(x,c,alpha)$lower
  }
#obtain upper end of CI
f <- function(theta) 2*pnorm(theta - c) -1 - (1 - alpha)* Q(theta)
theta2 <- uniroot( f,c( theta1,c + qnorm(1 - alpha/2) ) ) $ root      # theta2 - c cannot be greater than qnorm(1 - alpha/2)!
f <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) * Q(theta)
upper <- uniroot( f, c( theta2, x + 1.1 * qnorm(1 - alpha/2) ) ) $ root # upper is < x + 1.96, but for root-finding issues I put 2 * 1.96
CI <- list(lower = lower, upper = upper)
if (is.neg == 1) CI <- list(lower = -upper, upper = -lower)
return(CI)
}





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





QB.CI <- function(x,c,alpha){

f <- function(theta) 1 - pnorm(x - theta) - alpha / 2 * Q(theta)
lower <- uniroot(f, c(-c - qnorm(1 - alpha / 2), x - qnorm(1 - alpha / 2))) $ root
f <- function(theta) pnorm(theta - x) - (1 - alpha / 2) * Q(theta)
upper <- uniroot(f, c(x, x + qnorm(1 - alpha / 2))) $ root
CI <- list(lower = lower, upper = upper)
return(CI)
}





# Length of the QC CI at x = ctildezero.qc VS lambda

QC.maxlength <- function(lambda,c,alpha){

dmin <- function(theta) max( -c - theta + qnorm(1 - alpha * Q(theta)), 0 )
dmax <- function(theta) -c - Shortest.AR(theta,c,alpha)$A[1]
f <- function(d) 1 + lambda * ( 1 - dnorm(c + d) / dnorm(qnorm(2 - pnorm(c + d) - alpha * Q(0))) )
dminzero <- dmin(0)
dmaxzero <- dmax(0)
while (is.na(suppressWarnings(f(dminzero)))) dminzero <- dminzero + 1e-12  #just a technical issue with computing f() exactly at dminzero	
cbarzero.qc <- uniroot(f, c(dminzero, dmaxzero))$root + c

ctildezero.qc <- qnorm( 2 - pnorm(cbarzero.qc) - alpha * Q(0) )
a <- QC.CI(ctildezero.qc,lambda,c,alpha) $ lower
b <- QC.CI(ctildezero.qc,lambda,c,alpha) $ upper
l <- b - a
return(l)
}