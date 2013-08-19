## Post selection confidence intervals from:
# Benjamini, Y., and D. Yekutieli. 2005. 
# “False Discovery Rate-Adjusted Multiple Confidence Intervals for Selected Parameters.”
# Journal of the American Statistical Association 100 (469): 71–81.
# Construct an alpha level CI for a parameter of value x which has been selected as part
# of n.rejections in n.tests

BYCI<- function(x, sigsq, n.rejections, n.tests, alpha){
  alpha.adjusted<- 1 - n.rejections * alpha / (2*n.tests)
  arm<- sigsq * qnorm(alpha.adjusted)
  
  result<- c(ll=x-arm, ul=x+arm)
  return(result)    
}
## Testing:
# BYCI(x=1, sigsq=1, n.rejections=10, n.tests=100, alpha=0.05)
