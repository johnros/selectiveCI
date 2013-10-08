
pearson2zscore<- function(r, n){
  result<- NA
  if(!is.na(r)){
    stopifnot(abs(r) <= 1 && n%%1==0L && n>3 )
    
    rho.sd<- 1 / sqrt(n-3)
    result<- atanh(r)/rho.sd
  }
  return(result)
}
## Testing
#pearson2zscore(0.5, 400)
#pearson2zscore(NA, 400)



zscore2pearson<- function(z, n){
  result<- NA
  if(!is.na(z)){
    stopifnot( is.numeric(z) && n%%1==0L && n>3 )
    
    rho.sd<- 1 / sqrt(n-3)
    result<- tanh(z * rho.sd)
  }
  return(result)
}
## Testing:
#zscore2pearson(pearson2zscore(0.5, 400), 400)
#zscore2pearson(NA, 400)



## Check validity of inputs and return informative errors:
checkX<- function(x, cutoff){
  stopifnot(is.numeric(x))
  if(!missing(cutoff)){
    stopifnot(is.numeric(cutoff))
    if(abs(x) <= cutoff) stop('Observed value needs to be strictly greater than cutoff to be selected.')
  }
}

checkSigsq<- function(sigsq){
  stopifnot(is.numeric(sigsq))
  if(sigsq<0) stop('Variance cannot be negative. Check sigsq')
}

checkAlpha<- function(alpha){
  stopifnot(is.numeric(alpha))
  if(alpha<0 || alpha > 1) stop('Probabilities are between 0 and 1. Check alpha')
}

checkRatio<- function(ratio){
  stopifnot(is.numeric(ratio))
  if(ratio<0) stop('Invalid Modified Pratt ratio. Check r')
}

checkLambda<- function(lambda){
  stopifnot(is.numeric(lambda))
  if(lambda<0) stop('Negative penalties not allowed. Check lambda')
}

checkCutoff<- function(cutoff){
  stopifnot(is.numeric(cutoff))
  if(cutoff<0) stop('Selection cutoff needed in absolute value. Check cutoff')
}

checkTheta<- function(theta){
  stopifnot(is.numeric(theta))
}

