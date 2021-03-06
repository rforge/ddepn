# TODO: Add comment
# 
# Author: benderc
###############################################################################


#reverse.direction <- function(phi, i) {
#	cs <- coord(i,phi)
#	phi[cs[2],cs[1]] <- phi[cs[1],cs[2]]
#	phi[cs[1],cs[2]] <- 0
#	phi
#}

mcmc_move <- function(bestmodel, type) {
	stimuli <- bestmodel$stimuli
	gammaposs <- bestmodel$gammaposs
	nummoves <- bestmodel$nummoves
	numposs <- 0
	
	if(type=="addactivation" || type=="addinhibition") {
		phi <- bestmodel$phi
		diag(phi) <- 1
		phi[,unlist(stimuli)] <- matrix(1,nrow=nrow(phi),ncol=length(unlist(stimuli)))
		# exclude edges from nodes that are disconnected and thus 0 all the time
		noreach <- which(apply(gammaposs, 1, function(x) all(x==0)))
		for(i in noreach) {
			phi[i,] <- rep(1,ncol(phi))
		}
		poss <- which(phi==0)
		## pr that exactly the edge added is deleted, this is for the proposal density
		possback <- which(bestmodel$phi!=0)
		rm(phi)
	}
	if(type=="delete") {
		poss <- which(bestmodel$phi!=0)
		## pr that exactly the edge removed is added again, this is for the proposal density
		phi <- bestmodel$phi
		diag(phi) <- 1
		phi[,unlist(stimuli)] <- matrix(1,nrow=nrow(phi),ncol=length(unlist(stimuli)))
		# exclude edges from nodes that are disconnected and thus 0 all the time
		noreach <- which(apply(gammaposs, 1, function(x) all(x==0)))
		for(i in noreach) {
			phi[i,] <- rep(1,ncol(phi))
		}
		possback <- which(phi==0)
		
	}
	if(type=="switchtype") {
		poss <- which(bestmodel$phi!=0)
		possback <- poss
	}
	if(type=="revert") {
		poss <- which(bestmodel$phi!=0)
		for(p in poss) {
			cs <- coord(p,bestmodel$phi)
			if(bestmodel$phi[cs[2],cs[1]]!=0) {
				poss <- poss[-p]
			}
		}
		possback <- poss
	}
	# pr that edge is chosen in the selected move: P(Edge|move)
	pegm <- log2(1/(length(poss)*nummoves)) / max(length(poss),1)
	#pegm <- -log2(length(poss))
	# pr that the move from above for the given edge is reverted: P(Edgechangeundo|move)
	pegmundo <- log2(1/((length(possback)+1) * nummoves))  / max(length(possback),1)
	#pegmundo <- -log2((length(possback)+1))
	if(length(poss)>0) {
		tps <- bestmodel$tps
		stimuli <- bestmodel$stimuli
		reps <- bestmodel$reps
		datx <- bestmodel$dat
		maxiter <- bestmodel$maxiter
		TSA <- bestmodel$TSA
		Tt <- bestmodel$Tt
		lambda <- bestmodel$lambda
		B <- bestmodel$B
		Z <- bestmodel$Z
		# permute poss to propose moves in different orders
		#poss <- sample(poss)
		counter <- 1
		numbettermodel <- 1
		bettermodels <- list()
		
		# select the move randomly
		i <- ifelse(length(poss)>1,sample(poss, 1),poss[1])
		phi.n <- bestmodel$phi
		cds <- coord(i, phi.n)
		if(cds[1]==cds[2]) {
			print("ERROR: chosen element of the diagonal for the move.")
			browser()
		}
		switch(type,
				addactivation=phi.n[i] <- 1,
				addinhibition=phi.n[i] <- 2,
				switchtype=phi.n[i] <- phi.n[i]%%2 + 1,
				delete=phi.n[i] <- 0,
				reverse=phi.n <- reverse.direction(phi.n,i))	
		L.res <- perform.hmmsearch(phi.n, bestmodel)	
		gammaposs.n <- L.res$gammaposs
		gamma.n <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
		theta.n <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
		L.n <- L.res$Likl
		if(any(L.res$Likl==-Inf))
			browser()
		bic.n <- L.res$bic
		aic.n <- L.res$aic
		if(is.null(lambda)) {
			posterior.n <- NULL
		} else {
			posterior.n <- posterior(phi.n, L.n, lambda, B, Z)			
		}
		bettermodels[[numbettermodel]] <- list(phi=phi.n,L=L.n,aic=aic.n,bic=bic.n,posterior=posterior.n,dat=datx,
				theta=theta.n, gamma=gamma.n, gammaposs=gammaposs.n, tps=tps, stimuli=stimuli,
				reps=reps, maxiter=maxiter, TSA=NULL, Tt=NULL, lastmove=type, coords=cds,
				lambda=lambda, B=B, Z=Z, pegm=pegm, pegmundo=pegmundo,nummoves=bestmodel$nummoves)	
		numbettermodel <- numbettermodel + 1
	} else {
		bettermodels <- list(bestmodel)
	}
	return(bettermodels)
}

