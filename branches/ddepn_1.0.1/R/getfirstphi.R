# Helper for netga.R
# 
# Author: benderc
###############################################################################

getfirstphi <- function(x, datx, stimuli, V, tps, reps, maxiter, lambda=NULL, B=NULL, Z=NULL) {
	cat("-")
	if(is.null(x)) {
		phireference <- matrix(0,nrow=nrow(datx),ncol=nrow(datx),dimnames=list(rownames(datx),rownames(datx)))
		phiasis <- FALSE
	} else {
		phireference <- x
		phiasis <- TRUE
	}			
	# for all i>2 sample random start networks
	ps <- samplephi(phireference, stimuli, V, tps, reps, datx, searchstatespace=TRUE, maxiter=maxiter,
			phiasis=phiasis, lambda=lambda, B=B, Z=Z)
	model <- list(phi=ps$phi.n, L=ps$Lnew, aic=ps$aicnew, bic=ps$bicnew, posterior=ps$posterior,
			dat=datx, theta=ps$theta, gamma=ps$gammax, gammaposs=ps$gammaposs,tps=tps,
			stimuli=stimuli,reps=reps,maxiter=maxiter, TSA=NULL, Tt=NULL, lastmove=NULL,
			coords=c(1,1), lambda=lambda, B=B, Z=Z)
	model
}

