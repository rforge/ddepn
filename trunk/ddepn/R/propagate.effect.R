# Effect propagation for single- and multi-Experiment data.
# 
# Author: benderc
###############################################################################


propagate.effect.set <- function(phi,tps,stimuli,reps=1) {
	gam <- NULL
	for(i in 1:length(stimuli)) {
		stimulus <- stimuli[[i]]
		tmp <- propagate.effect.simple(phi,tps,stimulus,reps)		
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
propagate.effect.simple <- function(phi,tps,stimulus,reps=1) {
	gamma <- matrix(0,nrow=nrow(phi), ncol=length(tps), dimnames=list(rownames(phi),tps))
	# initialize the stimulated nodes
	gamma[stimulus,] <- 1
	stopprop <- FALSE
	tp <- 2
	states <- NULL
	while(!stopprop && tp!=ncol(gamma)) {
		for(prot in rownames(gamma)) {
			if(any(which(rownames(gamma)==prot) %in% stimulus))
				next
			status <- max(max(gamma[,(tp-1)] * phi[,prot]))
			# if inhibition found as input for the node, it will become inactive
			if(status==2) {
				status <- 0
			}
			gamma[prot,tp] <- status
		}
		states <- c(states,bin2dec(gamma[,tp]))
		tp <- tp + 1
		if(length(unique(states))==length(states[-length(states)])) {
			stopprop <- TRUE
			gamma <- gamma[,1:(tp-2)]
		}
	}
	if(reps>1) {
		gamma <- expand.gamma(gamma,reps)
	}
	gamma
}

