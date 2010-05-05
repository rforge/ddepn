# Simulate a random signalling network:
#   For a given number of stimuli, sample this number of stimuli nodes randomly.
#	Start at the stimuli and keep on adding activation edges to other random nodes
#	in V until all nodes were reached. Then add prop.inh*|E| random inhibition edges.
#   
#   n: number of nodes
#	nstim: number of stimuli
#	cstim: number of combinatorial stimuli, creates subsequently all pairs, triples, 
#          quadruples of stimuli combinations until number cstim is reached
#	prop.inh: proportion of inhibiting edges w.r.t. the number of activating edges
# Author: benderc
###############################################################################

signalnetwork_old <- function(n=10, nstim=2, cstim=0, prop.inh=.2,plot=F,gamma=1) {
	#V <- LETTERS[1:n]
	## initialise
	V <- paste("X",1:n,sep="")
	phi <- matrix(0,nrow=n, ncol=n, dimnames=list(V,V))
	stim <- sample(1:n, nstim)
	names(stim) <- V[stim]
	stimuli <- list()
	for(i in 1:length(stim)) {
		stimuli[[i]] <- stim[i]
	}
	# now the combinatorial stimuli
	if(cstim > 0) {
		lstim <- length(stimuli)
		k <- 2
		combs <- choose(nstim,k)
		while(combs < cstim) {
			k <- k+1
			if(k>nstim) {
				stop(paste("You want too many combinatorial stimuli. For ",nstim, "stimuli, only ", combs, " combinations are possible. Aborting."))
			}
			combs <- combs + choose(nstim,k)	
		}
		tuples <- singlestim <- names(unlist(stimuli))
		for(kprime in 2:k) {
			kprime <- kprime + 1
			newtuples <- NULL
			for(i in 1:nstim) {
				for(j in 1:length(tuples)) {
					if(!is.na(match(singlestim[i],strsplit(tuples[j],"-")[[1]])))
						next
					tup <- paste(singlestim[i], tuples[j], sep="-")
					newtuples <- c(newtuples, tup)
				}
			}
			newtuples <- unique(t(sapply(newtuples, function(x) paste(sort(strsplit(x,"-")[[1]]),sep="-"))))
			for(l in 1:nrow(newtuples)) {
				vec <- as.numeric(match(newtuples[l,],V))
				names(vec) <- as.character(newtuples[l,])
				stimuli[[(lstim+l)]] <- vec
				if(length(stimuli)==(nstim+cstim)) {
					kprime <- k+1
					break				
				}
			}
			lstim <- length(stimuli)
			tuples <- apply(newtuples, 1, paste, collapse="")
		}
	}
	
	# simulate the network now
	stimtmp <- unique(unlist(stimuli))
	ccl <- -1
	counter <- 1
	broke <- FALSE
	while(ccl!=1 && counter<=200) {
		cat(":")
		newstimtmp <- NULL
		for(st in stimtmp) {
			tosample <- setdiff(1:n,union(unlist(stimuli),stimtmp))
			if(length(tosample)<1) {
				degprob <- 0.432043005
				broke <- TRUE
				break
			} else {
				degprob <- (0:(length(tosample)-1))^(-gamma)
				degprob[1] = 1
				degprob = degprob/sum(degprob)
			}
			while(1==1) {
				cat(";")
				newedges <- NULL
				for(dp in 1:length(tosample)) {
					et <- sample(c(0,1), 1, prob=c((1-degprob[dp]), degprob[dp]))
					#if(et==1)
					#	et <- sample(c(1,2),1,prob=c((1-prop.inh),prop.inh)) # might be inhibition
					newedges <- c(newedges, et)
				}
					
					
				#newedges <- sample(c(0,1), length(tosample), prob=c(0.8,0.2), replace=T)
				if(!all(newedges==0) || length(newedges)==0)
					break
			}
			if(length(newedges)==0) {
				ccl<-1
				break
			}
			
			newstimtmp <- c(newstimtmp, tosample[newedges!=0])
			phi[st, tosample] <- newedges
			g <- as(phi, "graphNEL")
			ccl <- length(connComp(g))
			if(ccl==1)
				break
		}
		if(broke)
			break
		stimtmp <- unique(newstimtmp)
		counter <- counter + 1
	}
	phitmp <- phi
	diag(phitmp) <- 1
	phitmp[unique(unlist(stimuli)),] <- rep(1,ncol(phitmp))
	phitmp[,unique(unlist(stimuli))] <- rep(1,ncol(phitmp))
	remaining <- which(phitmp==0)
	theedges <- which(phi==1)
	numinh <- floor(length(theedges) * prop.inh)
	inhibs <- sample(remaining, numinh)
	phi[inhibs] <- 2
	if(plot)
		plotdetailed(phi,stimuli=stimuli)
	return(list(phi=phi, stimuli=stimuli))
}


signalnetwork <- function(n=10, nstim=2, cstim=0, prop.inh=.2,plot=F,gamma=1,B=NULL) {
	#V <- LETTERS[1:n]
	## initialise
	V <- paste("X",1:n,sep="")
	phi <- matrix(0,nrow=n, ncol=n, dimnames=list(V,V))
	stim <- sample(1:n, nstim)
	names(stim) <- V[stim]
	stimuli <- list()
	for(i in 1:length(stim)) {
		stimuli[[i]] <- stim[i]
	}
	# now the combinatorial stimuli
	if(cstim > 0) {
		lstim <- length(stimuli)
		k <- 2
		combs <- choose(nstim,k)
		while(combs < cstim) {
			k <- k+1
			if(k>nstim) {
				stop(paste("You want too many combinatorial stimuli. For ",nstim, "stimuli, only ", combs, " combinations are possible. Aborting."))
			}
			combs <- combs + choose(nstim,k)	
		}
		tuples <- singlestim <- names(unlist(stimuli))
		for(kprime in 2:k) {
			kprime <- kprime + 1
			newtuples <- NULL
			for(i in 1:nstim) {
				for(j in 1:length(tuples)) {
					if(!is.na(match(singlestim[i],strsplit(tuples[j],"-")[[1]])))
						next
					tup <- paste(singlestim[i], tuples[j], sep="-")
					newtuples <- c(newtuples, tup)
				}
			}
			newtuples <- unique(t(sapply(newtuples, function(x) paste(sort(strsplit(x,"-")[[1]]),sep="-"))))
			for(l in 1:nrow(newtuples)) {
				vec <- as.numeric(match(newtuples[l,],V))
				names(vec) <- as.character(newtuples[l,])
				stimuli[[(lstim+l)]] <- vec
				if(length(stimuli)==(nstim+cstim)) {
					kprime <- k+1
					break				
				}
			}
			lstim <- length(stimuli)
			tuples <- apply(newtuples, 1, paste, collapse="")
		}
	}
#browser()
	# simulate the network now
	stimtmp <- unique(unlist(stimuli))
	ccl <- -1
	counter <- 1
	broke <- FALSE
	while(ccl!=1 && counter<=200) {
		cat(":")
		newstimtmp <- NULL
		for(st in stimtmp) {
			tosample <- setdiff(1:n,union(unlist(stimuli),stimtmp))
			if(length(tosample)<1) {
				degprob <- 0.432043005
				broke <- TRUE
				break
			} else {
				degprob <- (0:(length(tosample)-1))^(-gamma)
				degprob[1] = 1
				degprob = degprob/sum(degprob)
			}
			while(1==1) {
				cat(";")
				newedges <- NULL
				for(dp in 1:length(tosample)) {
					et <- sample(c(0,1), 1, prob=c((1-degprob[dp]), degprob[dp]))
					#if(et==1)
					#	et <- sample(c(1,2),1,prob=c((1-prop.inh),prop.inh)) # might be inhibition
					newedges <- c(newedges, et)
				}
				
				
				#newedges <- sample(c(0,1), length(tosample), prob=c(0.8,0.2), replace=T)
				if(!all(newedges==0) || length(newedges)==0)
					break
			}
			if(length(newedges)==0) {
				ccl<-1
				break
			}
			if(!is.null(B)) {
				# now it is clear that a certain number of edges is sampled
				# i.e. the number of 1s in newedges
				# now find out to which node the edge is going to using 
				# the probabilities of the B-prior-matrix
				probs <- (B[st,]/sum(B[st,]))[tosample]
				nr <- length(which(newedges!=0))
				edd <- sample(tosample, nr, prob=probs)
			} else {
				edd <- tosample[newedges!=0]
			}
			newstimtmp <- c(newstimtmp, edd)
			phi[st, tosample] <- newedges
			g <- as(phi, "graphNEL")
			ccl <- length(connComp(g))
			if(ccl==1)
				break
		}
		if(broke)
			break
		stimtmp <- unique(newstimtmp)
		counter <- counter + 1
	}
	phitmp <- phi
	diag(phitmp) <- 1
	phitmp[unique(unlist(stimuli)),] <- rep(1,ncol(phitmp))
	phitmp[,unique(unlist(stimuli))] <- rep(1,ncol(phitmp))
	remaining <- which(phitmp==0)
	theedges <- which(phi==1)
	numinh <- floor(length(theedges) * prop.inh)
	inhibs <- sample(remaining, numinh)
	phi[inhibs] <- 2
	if(plot)
		plotdetailed(phi,stimuli=stimuli)
	return(list(phi=phi, stimuli=stimuli))
}
