# AIC and BIC scores.
# 
# Author: benderc
###############################################################################


get.aic <- function(phi,L) {
	k <- length(which(phi!=0))
	aic <- 2*k - 2*L
	aic
}
get.bic <- function(phi,L,n) {
	k <- length(which(phi!=0))
	bic <- -2*L + k*log2(n)
	bic
}

