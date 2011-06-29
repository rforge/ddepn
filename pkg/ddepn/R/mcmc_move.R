# TODO: Add comment
# 
# Author: benderc
###############################################################################

mcmc_move <- function(bestmodel, type) {
	stimuli <- bestmodel$stimuli
	gammaposs <- bestmodel$gammaposs
	nummoves <- bestmodel$nummoves
	phiorig <- bestmodel$phi.orig
	scale_lik <- bestmodel$scale_lik
	numposs <- 0
	fanin <- bestmodel$fanin
	if(type %in% c("addinhibition","addactivation","add")) {
		phi <- bestmodel$phi
		diag(phi) <- 1
		phi[,unlist(stimuli)] <- matrix(1,nrow=nrow(phi),ncol=length(unlist(stimuli)))
		# exclude edges from nodes that are disconnected and thus 0 all the time
		# as well as edges to nodes having already more than 4 incoming edges
		fanin_omit <- which(colSums(detailed.to.simple.regulations(bestmodel$phi))>=fanin)
		for(i in fanin_omit) {
			phi[,i] <- rep(1,nrow(phi))
		}
		noreach <- which(apply(gammaposs, 1, function(x) all(x==0)))
		for(i in noreach) {
			phi[i,] <- rep(1,ncol(phi))
		}
		poss <- which(phi==0)
		## pr that exactly the edge added is deleted, this is for the proposal density
		possback <- which(bestmodel$phi!=0) + 1 # all edges that are not 0 plus the one added
		rm(phi)
#		if(length(poss)>0) {
#			## check if edges would point to a stimulus node
#			cds <- sapply(poss, coord, mat=bestmodel$phi)
#			if(any(cds[2,] %in% unique(unlist(stimuli)))) {
#				print("Error: add edge would be drawn to a stimulus.")
#				browser()
#			}
#		}
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
	if(type=="revert" | type=="revswitch") {
		phi <- bestmodel$phi
		## which of the nodes already have fanin incoming edges
		fanin_omit <- which(colSums(detailed.to.simple.regulations(bestmodel$phi))>=fanin)
		## set all outgoing edges of these nodes to 0, since they must not be reverted
		for(i in fanin_omit) {
			phi[i,] <- rep(0,ncol(phi))
		}
		## all edges that would go into a stimulus if reverting are not allowed to revert
		phi[unique(unlist(stimuli)), ] <- 0
		## exclude edges that are disconnected
		#noreach <- which(apply(gammaposs, 1, function(x) all(x==0)))
		#for(i in noreach) {
		#	phi[i,] <- rep(1,ncol(phi))
		#}
		
		## which edges are allowed to revert/revswitch
		poss <- which(phi!=0)
		for(p in poss) {
			cs <- coord(p,bestmodel$phi)
			## do not allow to revert/revswitch if this edge already exists
			if(bestmodel$phi[cs[2],cs[1]]!=0) {
				poss <- poss[-p]
			}
		}
		possback <- poss
#		if(length(poss)>0) {			
#			## check if edges would point to a stimulus node
#			cds <- sapply(poss, coord, mat=bestmodel$phi)
#			if(any(cds[1,] %in% unique(unlist(stimuli)))) {
#				print("Error: revswitch edge would be drawn to a stimulus.")
#				browser()
#			}
#		}
	}
	if(length(poss)>0) {		
		## pr that edge is chosen in the selected move: P(Edge|move)
		pegm <- -log(length(poss)) - log(nummoves) ## v3
		## pr that the move from above for the given edge is reverted: P(Edgechangeundo|move)
		pegmundo <- -log(length(possback)) - log(nummoves) ## v3
		## get the model variables
		tps <- bestmodel$tps
		stimuli <- bestmodel$stimuli
		reps <- bestmodel$reps
		dat <- bestmodel$dat
		hmmiterations <- bestmodel$hmmiterations
		TSA <- bestmodel$TSA
		Tt <- bestmodel$Tt
		lambda <- bestmodel$lambda
		B <- bestmodel$B
		Z <- bestmodel$Z
		gam <- bestmodel$gam
		it <- bestmodel$it
		K <- bestmodel$K
		priortype <- bestmodel$priortype
		allow.stim.off <- bestmodel$allow.stim.off
		counter <- 1
		numbettermodel <- 1
		bettermodels <- list()
		
		## select the move randomly
		i <- ifelse(length(poss)>1,sample(poss, 1),poss[1])
		phi.n <- bestmodel$phi
		cds <- coord(i, phi.n)
		if(cds[1]==cds[2] || cds[2] %in% unique(unlist(stimuli))) {
			print("ERROR: chosen element of the diagonal for the move or edge back to stimulus.")
			browser()
		}
		switch(type,
				add=phi.n[i] <- sample(c(1,2),1),
				addactivation=phi.n[i] <- 1,
				addinhibition=phi.n[i] <- 2,
				switchtype=phi.n[i] <- phi.n[i]%%2 + 1,
				delete=phi.n[i] <- 0,
				reverse=phi.n <- reverse.direction(phi.n,i),
				revswitch=phi.n <- reverse.direction(phi.n,i,switchtype=TRUE))
		## debug me if constraints are violated
		if(any(colSums(detailed.to.simple.regulations(phi.n))>fanin)) {
			print("Oops, mcmc_move seems to produce more incoming edges than allowed, according to fanin-setting.")
			browser()
		}
		## HMM for the proposed network
		L.res <- perform.hmmsearch(phi.n, bestmodel)	
		gammaposs.n <- L.res$gammaposs
		gamma.n <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
		theta.n <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
		L.n <- L.res$Likl
		if(any(L.res$Likl==-Inf))
			browser()
		bic.n <- L.res$bic
		aic.n <- L.res$aic
		pr.n <- prior(phi.n, lambda, B, Z, gam, it, K, priortype)
		if(priortype %in% c("laplaceinhib", "laplace", "scalefree", "uniform")) {
			posterior.n <- L.n + pr.n
		} else {
			posterior.n <- NULL			
		}
		bettermodels[[numbettermodel]] <- list(phi=phi.n,L=L.n,aic=aic.n,bic=bic.n,posterior=posterior.n,dat=dat,
				theta=theta.n, gamma=gamma.n, gammaposs=gammaposs.n, tps=tps, stimuli=stimuli,
				reps=reps, hmmiterations=hmmiterations, TSA=NULL, Tt=NULL, lastmove=type, coords=cds,
				lambda=lambda, B=B, Z=Z, pegm=pegm, pegmundo=pegmundo,nummoves=bestmodel$nummoves,fanin=fanin,
				gam=gam, it=it, K=K,phi.orig=phiorig,priortype=priortype,pr=pr.n,
				mu_run=bestmodel$mu_run,Qi=bestmodel$Qi,sd_run=bestmodel$sd_run, scale_lik=scale_lik, allow.stim.off=allow.stim.off)
				#,mean_thetax=bestmodel$mean_thetax, mean_squared_thetax=bestmodel$mean_squared_thetax,
				#sd_thetax=bestmodel$sd_thetax)	
		numbettermodel <- numbettermodel + 1
	} else {
		bettermodels <- NULL
	}
	return(bettermodels)
}
