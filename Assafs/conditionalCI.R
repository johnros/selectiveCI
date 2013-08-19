#### Compute conditional CI ####

## Sketch:
# Find distribution of data conditional on rejection
# Find highest density CIs 

conditional.pnorm<- function(x, cutoff, mean=0, sd=1){
	## Initializing:
	stopifnot(cutoff>0)
	result<- NA
	scaler<- pnorm(abs(cutoff), mean=mean, sd=sd)-pnorm(-abs(cutoff), mean=mean, sd=sd)

	if(x < cutoff && abs(x) > cutoff) {
		result<- pnorm(x, mean=mean, sd=sd)/(1-scaler)
	}

	else if(x > cutoff) {
		result<- (pnorm(x, mean=mean, sd=sd) - scaler) /(1-scaler)
	}

	else if(abs(x) <= cutoff) {
		result<- pnorm(-abs(cutoff), mean=mean, sd=sd)/(1-scaler)
	}
	return(result)
}
## Testing:
conditional.pnorm(x=-1, cutoff=0.52, mean=0, sd=1)
conditional.pnorm(x=-0.2, cutoff=0.2, mean=1, sd=1)
conditional.pnorm(x=0.2, cutoff=0.2, mean=1, sd=1)
conditional.pnorm(x=-2, cutoff=2, mean=0, sd=1)
x<- seq(-5,5, length=1000)
plot(Vectorize(conditional.pnorm)(x=x, cutoff=1, mean=1, sd=1)~ x, cex=0.1)


## Find the quantiles of the conditional Gaussian:
conditional.qnorm<- function(q, cutoff, mean=0, sd=1){
	stopifnot(cutoff>0 && q<= 1 && q>= 0)

	if(q==1) return(Inf)
	if(q==0) return(-Inf)

	cutoff.probability<- conditional.pnorm(x=cutoff, cutoff=cutoff, mean=mean, sd=sd)
	result<- NA
	q.lower<- ifelse(q < 0.5 , q, 1-q)

	# Target function:
	.function<- function(x) conditional.pnorm(x, cutoff=cutoff, mean=mean, sd=sd)-q

	if(q < cutoff.probability) .interval <- c(qnorm(q.lower, mean, 5*sd), -cutoff)
	else .interval <- c(cutoff, qnorm(1-q.lower, mean, 5*sd))
	try(result<- uniroot(f=.function, 
			     interval=.interval )$root)     
	return(result)
}
## Testing:
conditional.qnorm(0.2, 1)
conditional.qnorm(0.4, 1)
conditional.qnorm(0.001, 1)



## HDR ##
conditional.HDR<- function(alpha, cutoff, mean, sd){
  stopifnot(alpha<0.5 && cutoff >0)
  # alpha<- 0.1
  # cutoff<- 1
  # mean<- 1
  # sd<- 1
	center.dist<- abs(mean) - cutoff
  result<- list(L=NA, H=NA)
  
	A.probability<- 
    conditional.pnorm(mean+center.dist, cutoff=cutoff, mean=mean, sd=sd) 
	  - conditional.pnorm(mean-center.dist, cutoff=cutoff, mean=mean, sd=sd)

  if(1-alpha < A.probability){
    lower.quantile<- conditional.qnorm(q=alpha, cutoff=cutoff, mean=mean, sd=sd)
    upper.quantile<- mean+(mean-lower.quantile)
    result$L<- min(lower.quantile, upper.quantile)
    result$H<- max(lower.quantile, upper.quantile)
  }
  
  else{
    center.dist2<- mean + cutoff
    scaler<- pnorm(abs(cutoff), mean=mean, sd=sd)-pnorm(-abs(cutoff), mean=mean, sd=sd)
    B.probability<- scaler/(1-scaler)    
    if(1-alpha < A.probability + B.probability){
      lower.quantile<- cutoff
      upper.quantile<- 
        conditional.qnorm(q= conditional.pnorm(x=-cutoff, cutoff=cutoff, mean=mean, sd=sd) 
                          +(1-alpha), cutoff=cutoff, mean=mean, sd=sd)
      result$L<- lower.quantile
      result$H<- upper.quantile
    }
    else{
      NULL
    }
  }
      
  return(result)  
}
## Testing:
conditional.HDR(alpha=0.1, cutoff=1, mean=1, sd=1)





## HDR CI ##:
conditional.CI<- function(x, cutoff, sd){

}

## Equi-probability CI ##
