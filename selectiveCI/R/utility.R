
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