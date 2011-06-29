# Effect propagation for single- and multi-Experiment data.
# 
# Author: benderc
###############################################################################


propagate.effect.set <- function(phi,stimuli,allow.stim.off=FALSE) {
	gam <- NULL
	for(i in 1:length(stimuli)) {
		stimulus <- stimuli[[i]]
		tmp <- propagate.effect.simple(phi,stimulus,stimuli,allow.stim.off=allow.stim.off)		
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
##
## now also allows for the stimulus to loose activity, set by argument allow.stim.off
propagate.effect.simple <- function (phi, stimulus, stimuli, allow.stim.off=FALSE) 
{
	gamma <- matrix(0, nrow = nrow(phi), ncol = 1, dimnames = list(rownames(phi), 
					1))
	gamma[stimulus, ] <- 1
	stopprop <- FALSE
	tp <- 1
	states <- bin2dec(gamma[, 1])
	while (!stopprop) {
		newcolumn <- rep(0, nrow(phi))
		newcolumn[stimulus] <- 1
		for (k in 1:nrow(gamma)) {
			prot <- rownames(gamma)[k]
			if (any(which(rownames(gamma) == prot) %in% unique(unlist(stimuli)))) 
				next
			status <- max(gamma[, tp] * phi[, prot])
			if (status == 2) {
				status <- 0
			}
			newcolumn[k] <- status
		}
		states <- c(states, bin2dec(newcolumn))
		tp <- tp + 1
		if (length(unique(states)) == length(states[-length(states)])) {
			stopprop <- TRUE
		}
		else {
			gamma <- cbind(gamma, newcolumn)
			colnames(gamma) <- 1:ncol(gamma)
		}
	}
	if(allow.stim.off) {
		M <- ncol(gamma)
		gamma <- cbind(gamma,gamma)
		colnames(gamma) <- 1:ncol(gamma)
		gamma[stimulus,(M+1):ncol(gamma)] <- 0
	}
	gamma
}
