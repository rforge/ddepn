# Helper for netga.R
# 
# Author: benderc
###############################################################################

getfirstphi <- function(x,dat,stimuli,V,tps,reps,hmmiterations,
		lambda=NULL,B=NULL,Z=NULL,fanin=4,gam=NULL,it=NULL,K=NULL) {
	cat("-")
	if(is.null(x)) {
		phireference <- matrix(0,nrow=nrow(dat),ncol=nrow(dat),dimnames=list(rownames(dat),rownames(dat)))
		phiasis <- FALSE
	} else {
		phireference <- x
		phiasis <- TRUE
	}			
	# for all i>2 sample random start networks
	ps <- samplephi(phireference,stimuli,V,tps,reps,dat,searchstatespace=TRUE,hmmiterations=hmmiterations,
			phiasis=phiasis,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K)
	model <- list(phi=ps$phi.n, L=ps$Lnew, aic=ps$aicnew, bic=ps$bicnew, posterior=ps$posterior,
			dat=dat, theta=ps$theta, gamma=ps$gammax, gammaposs=ps$gammaposs,tps=tps,
			stimuli=stimuli,reps=reps,hmmiterations=hmmiterations, TSA=NULL, Tt=NULL, lastmove=NULL,
			coords=c(1,1), lambda=lambda, B=B, Z=Z, fanin=fanin, gam=gam, it=it, K=K)
	model
}

