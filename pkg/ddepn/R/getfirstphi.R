# Helper for netga.R
# 
# Author: benderc
###############################################################################

getfirstphi <- function(x,dat,stimuli,V,tps,reps,hmmiterations,
		lambda=NULL,B=NULL,Z=NULL,fanin=4,gam=NULL,it=NULL,K=NULL,priortype="none",scale_lik=FALSE,allow.stim.off=TRUE,
		implementation="C") {
	cat("-")
	if(is.null(x)) {
		phireference <- matrix(0,nrow=nrow(dat),ncol=nrow(dat),dimnames=list(rownames(dat),rownames(dat)))
		phiasis <- FALSE
	} else {
		phireference <- x
		phiasis <- TRUE
	}			
	# get a random network
	ps <- samplephi(phireference,stimuli,V,tps,reps,dat,searchstatespace=TRUE,hmmiterations=hmmiterations,
			phiasis=phiasis,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K,priortype=priortype, scale_lik=scale_lik,allow.stim.off=allow.stim.off,
			implementation=implementation)
	## now create the initial model object, holding information on a specific network
	model <- list(phi=ps$phi.n, L=ps$Lnew, aic=ps$aicnew, bic=ps$bicnew, posterior=ps$posterior,
			dat=dat, theta=ps$theta, gamma=ps$gammax, gammaposs=ps$gammaposs,tps=tps,
			stimuli=stimuli,reps=reps,hmmiterations=hmmiterations, lastmove=NULL,
			coords=c(1,1), lambda=lambda, B=B, Z=Z, fanin=fanin, gam=gam, it=it, K=K, pr=ps$pr, priortype=priortype,scale_lik=scale_lik,
			allow.stim.off=allow.stim.off,implementation=implementation)
	model
}

