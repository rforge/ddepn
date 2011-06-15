# Effect propagation for single- and multi-Experiment data.
# 
# Author: benderc
###############################################################################


propagate.effect.set <- function(phi,stimuli) {
	gam <- NULL
	for(i in 1:length(stimuli)) {
		stimulus <- stimuli[[i]]
		tmp <- propagate.effect.simple(phi,stimulus,stimuli)		
		colnames(tmp) <- paste(paste(names(stimulus),collapse="&"), colnames(tmp), sep="_")
		gam <- cbind(gam, tmp)
	}
	gam
}

bin2dec <- function(x) {	
	sum(x * 2^(rev(seq(along=x)) - 1))
}

## propagate the input effects through a given network phi
## stop propagation as soon as a state is found twice
##
## note1: argument tps gives the length of the matrix to be calculated: experimenting with
## dynamically building gamma, so the number of columns (i.e. tps) doesn't have to be known
## before. hope this will not slow down too much if the number of nodes gets large... However,
## initializing a gamma matrix with 2^N columns is not possible for say N>20, which is bad.
## TODO: check if tps can be removed
##
## note2: argument reps doesn't seem to be necessary, since the propagation should create 
## unique states. TODO: check if this can be removed
propagate.effect.simple <- function(phi,stimulus,stimuli) {
	gamma <- matrix(0,nrow=nrow(phi), ncol=1, dimnames=list(rownames(phi),1))
	# initialize the stimulated nodes
	gamma[stimulus,] <- 1
	stopprop <- FALSE
	tp <- 1 ## if maxcolumns reached, stop
	states <- bin2dec(gamma[,1])
	while(!stopprop) {
		newcolumn <- rep(0,nrow(phi))
		newcolumn[stimulus] <- 1
		#for(prot in rownames(gamma)) {
		for(k in 1:nrow(gamma)) {
			prot <- rownames(gamma)[k]
			#if(any(which(rownames(gamma)==prot) %in% stimulus))
			if(any(which(rownames(gamma)==prot) %in% unique(unlist(stimuli))))
				next
			#status <- max(max(gamma[,(tp-1)] * phi[,prot])) # why double max?
			status <- max(gamma[,tp] * phi[,prot])
			# if inhibition found as input for the node, it will become inactive
			if(status==2) {
				status <- 0
			}
			newcolumn[k] <- status
			#newcolumn[tp] <- status
		}
		states <- c(states,bin2dec(newcolumn))
		tp <- tp + 1
		## if state was generated that is already in gamma, stop
		if(length(unique(states))==length(states[-length(states)])) {
			stopprop <- TRUE
		} else {
			## include the new state into gamma
			gamma <- cbind(gamma, newcolumn)
			colnames(gamma) <- 1:ncol(gamma)
		}
	}
	gamma
}
#propagate.effect.set.old <- function(phi,tps,stimuli,reps=rep(1,length(stimuli))) {
#	gam <- NULL
#	for(i in 1:length(stimuli)) {
#		stimulus <- stimuli[[i]]
#		tmp <- propagate.effect.simple(phi,tps,stimulus,reps[i])		
#		colnames(tmp) <- paste(paste(names(stimulus),collapse="&"), colnames(tmp), sep="_")
#		gam <- cbind(gam, tmp)
#	}
#	gam
#}
#propagate.effect.simple.old <- function(phi,tps,stimulus,reps=1) {
#	#gamma <- matrix(0,nrow=nrow(phi), ncol=length(tps), dimnames=list(rownames(phi),tps))
#	gamma <- matrix(0,nrow=nrow(phi), ncol=2, dimnames=list(rownames(phi),1:2))
#	# initialize the stimulated nodes
#	gamma[stimulus,] <- 1
#	stopprop <- FALSE
#	tp <- 2 ## if maxcolumns reached, stop
#	states <- NULL
#	while(!stopprop && tp!=ncol(gamma)) {
#		for(prot in rownames(gamma)) {
#			if(any(which(rownames(gamma)==prot) %in% stimulus))
#				next
#			#status <- max(max(gamma[,(tp-1)] * phi[,prot])) # why double max?
#			status <- max(gamma[,(tp-1)] * phi[,prot])
#			# if inhibition found as input for the node, it will become inactive
#			if(status==2) {
#				status <- 0
#			}
#			gamma[prot,tp] <- status
#		}
#		states <- c(states,bin2dec(gamma[,tp]))
#		tp <- tp + 1
#		if(length(unique(states))==length(states[-length(states)])) {
#			stopprop <- TRUE
#			gamma <- gamma[,1:(tp-2)]
#		}
#		newcolumn <- rep(0,nrow(phi))
#		newcolumn[stimulus] <- 1
#		gamma <- cbind(gamma, newcolumn)
#		colnames(gamma) <- 1:ncol(gamma)
#	}
#	if(reps>1) {
#		gamma <- expand.gamma(gamma,reps)
#	}
#	gamma
#}

