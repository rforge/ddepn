# Guess what
# 
# Author: benderc
###############################################################################


get.aic <- function(phi,L) {
	#diag(phi) <- 1
	#k <- length(which(phi==1))
	k <- length(which(phi!=0))
	#diag(phi) <- 0
	aic <- 2*k - 2*L
	aic
}
get.bic <- function(phi,L,n) {
	#diag(phi) <- 1
	#k <- length(which(phi==1))
	k <- length(which(phi!=0))
	#diag(phi) <- 0
	bic <- -2*L + k*log2(n)
	bic
}

