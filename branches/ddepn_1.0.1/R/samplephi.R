# TODO: Add comment
# 
# Author: benderc
###############################################################################


samplephi <- function(phi,stimuli, antibodies, tps, reps, dat, searchstatespace=FALSE,
		maxiter=5, phiasis=FALSE, lambda=NULL, B=B, Z=Z) {
	if(phiasis) {
		phi.n <- phi
	} else {
		#phi <- matrix(0, nrow=length(antibodies), ncol=length(antibodies), dimnames=list(antibodies,antibodies))
		nostimmat <- phi[,-unique(unlist(stimuli))]
		proplength <- length(antibodies)-length(unique(unlist(stimuli)))
		propnames <- colnames(nostimmat)
		prop <- matrix(sample(c(0,1,2),length(nostimmat),replace=TRUE),
				nrow=nrow(nostimmat),ncol=ncol(nostimmat), dimnames=dimnames(nostimmat))
		for(p in propnames)
			prop[p,p] <- 0
		phi.n <- cbind(phi[,unique(unlist(stimuli)),drop=F],prop)
		phi.n <- phi.n[,colnames(phi)]
	}
	longprop <- 1:max(length(tps),(nrow(phi.n)*100))
	gammaposs <- propagate.effect.set(phi.n,longprop,stimuli,reps=reps)
	gammaposs <- uniquegammaposs(gammaposs)
	# now get an initial gamma matrix
	gammax <- NULL
	for(sti in 1:length(stimuli)) {
		st <- stimuli[[sti]]
		indices <- grep(paste("^",names(st),"&",sep=""),colnames(gammaposs))
		gx <- replicatecolumns(gammaposs[,sort(sample(indices,length(tps),replace=TRUE))],reps)
		gammax <- cbind(gammax, gx)
	}
	#gammax <- replicatecolumns(gammaposs[,sort(sample(ncol(gammaposs),length(tps),replace=TRUE))],reps)
	#gammax <- propagate.effect.set(phi.n,tps,stimuli,reps=reps)
	Ltmplist <- likl(dat,gammax)
	Ltmp <- Ltmplist$L
	thetax <- Ltmplist$theta
	Ltmp[Ltmp==Inf] <- 0
	Ltmp[Ltmp==-Inf] <- 0
	Lnew <- sum(Ltmp)
	bicnew <- get.bic(phi.n,Lnew,length(dat))
	aicnew <- get.aic(phi.n,Lnew)
	if(searchstatespace) {
		bestmodel <- list(phi=phi.n,L=Lnew,aic=aicnew,bic=bicnew,dat=dat,
				theta=thetax, gamma=gammax, gammaposs=gammaposs, tps=tps, stimuli=stimuli, reps=reps,
				maxiter=maxiter, TSA=NULL, Tt=NULL, lastmove="addactivation", coords=c(1,1))
		L.res <- perform.hmmsearch(phi.n, bestmodel)	
		gammax <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
		thetax <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
		Lnew <- L.res$Likl
		bicnew <- L.res$bic
		aicnew <- L.res$aic		
	}
	# if prior is given
	if(is.null(lambda)) {
		postnew <- NULL
	} else {
		postnew <- posterior(phi.n, Lnew, lambda, B, Z)
	}
	return(list(phi.n=phi.n,gammax=gammax,thetax=thetax,posterior=postnew,Lnew=Lnew,bicnew=bicnew,aicnew=aicnew,gammaposs=gammaposs,lambda=lambda,B=B,Z=Z))
}


initialphi <- function(dat, phi, stimuli, Lmax, thetax, gammax, gammaposs,
		tps, reps, antibodies, n=100, multicores=FALSE, lambda=NULL, B=NULL, Z=NULL) {
	phimax <- phi
	thetamax <- thetax
	gammamax <- gammax
	gammaposs <- gammaposs
	bicmin <- get.bic(phi,Lmax,length(dat))
	aicmin <- get.aic(phi,Lmax)
	jobs <- list()
	for(i in 1:n) {
		if(i%%10==1)
			cat(".")
		
		if(multicores) {
			jobs[[i]] <- parallel(samplephi(phimax,stimuli, antibodies, tps, reps, dat, lambda=lambda, B=B, Z=Z))	
		} else {
			jobs[[i]] <- samplephi(phimax,stimuli, antibodies, tps, reps, dat, lambda=lambda, B=B, Z=Z)
		}
	}
	if(multicores)
		jobs.res <- collect(jobs, wait=TRUE)
	else
		jobs.res <- jobs
	
	for(i in 1:length(jobs.res)) {
		res <- jobs.res[[i]]
		if(res$bicnew<bicmin) {
			print(paste("Improved: ",res$bicnew))
			phimax <- res$phi.n
			Lmax <- res$Lnew
			aicmin <- res$aicnew
			bicmin <- res$bicnew
			posteriormax <- res$posterior
			thetamax <- res$thetax
			gammamax <- res$gammax
			gammaposs <- res$gammaposs
		}
	}
	return(list(phi=phimax, L=Lmax, aic=aicmin, bic=bicmin, posterior=posteriormax, thetax=thetamax, gammax=gammamax, gammaposs=gammaposs, lambda=lambda, B=B, Z=Z))
}

