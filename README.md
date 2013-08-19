Selective Confidence Intervals
===========

Compute confidence intervals for selected parameters as described in [1] and [2].




Reference
----------------------
[1] Weinstein, Asaf, William Fithian, and Yoav Benjamini. “Selection Adjusted Confidence Intervals With More Power to Determine the Sign.” Journal of the American Statistical Association 108, no. 501 (2013): 

[2] Benjamini, Y., and D. Yekutieli. 2005. “False Discovery Rate-Adjusted Multiple Confidence Intervals for Selected Parameters.” Journal of the American Statistical Association 100 (469): 71–81.



Installation
-------------
```{r}
library(devtools)
install_github("selectiveCI", "johnros", subdir='selectiveCI')
library(selectiveCI)
