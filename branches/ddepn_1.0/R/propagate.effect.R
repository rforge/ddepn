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

## propagate the input effects through a given network phi
## returns a matrix containing the stati of each proteins at each 
## timepoint
## all steps equally fast, once activated and not inhibited, a node
## stays active forever
propagate.effect.simple <- function(phi,tps,stimulus,reps=1) {
	gamma <- matrix(0,nrow=nrow(phi), ncol=length(tps), dimnames=list(rownames(phi),tps))
	# initialize the stimulated nodes
	#gamma[stimulus,1] <- 1
	gamma[stimulus,] <- 1
	
	# now propagate the effect through the network
	for(tp in 2:ncol(gamma)) {
		for(prot in rownames(gamma)) {
			#browser()
			if(any(which(rownames(gamma)==prot) %in% stimulus))
				next
			#status <- max(gamma[prot,(tp-1)],max(gamma[,(tp-1)] * phi[,prot]))
			status <- max(max(gamma[,(tp-1)] * phi[,prot]))
			# if inhibition found as input for the node, it will become inactive
			if(status==2) {
				status <- 0
			}
			gamma[prot,tp] <- status
		}
	}
	if(reps>1) {
		gamma <- expand.gamma(gamma,reps)
	}
	gamma
}


### propagate the input effects through a given network phi
### returns a matrix containing the stati of each proteins at each 
### timepoint
### activation happens immediately, inhibition has a delay of 1 step,
### activated proteins become inactive after a delay of 2 timesteps 
### after lost activation
#propagate.effect.simple2 <- function(phi.orig,tps,stimulus,reps=1) {
#	maxtp <- length(tps)
#	phi <- cbind(rep(0,nrow(phi.orig)),phi.orig)
#	phi <- rbind(rep(0,ncol(phi)),phi)
#	rownames(phi) <- colnames(phi) <- c("S",rownames(phi.orig))
#	phi[1,stimulus+1] <- 1
#	
#	gamma <- matrix(0,nrow=nrow(phi), ncol=maxtp+1, dimnames=list(rownames(phi),1:(maxtp+1)))
#	gamma["S",] <- 1
#	gamma[stimulus+1,1] <- 1
#	tp <- 1
#	msg <- matrix(0,nrow=nrow(phi),ncol=ncol(phi),dimnames=dimnames(phi))
#	# activations
#	tmp <- which((phi*gamma[,tp])==1)
#	targets.coords <- cbind((tmp%%nrow(phi)),(floor(tmp/nrow(phi)))+1)
#	## which ones must be set
#	msg[targets.coords[which(msg[targets.coords]==0),]] <- 1
#	
#	# inhibitions
#	tmp <- which((phi*gamma[,tp])==2)
#	targets.coords <- cbind((tmp%%nrow(phi)),(floor(tmp/nrow(phi)))+1)
#	## which ones must be set
#	msg[targets.coords[which(msg[targets.coords]==0),]] <- 2		
#	msg.new <- msg
#	
#	tp <- 2
#	while(tp<maxtp+1) {
#		msg <- msg.new
#		# decrease all messages>0 by 1, i.e. all effects connected to
#		# edges that become active some time are handled
#		msg.new[which(msg.new>0)] <- msg.new[which(msg.new>0)]-1
#		# increase all messages<0 by 1, i.e. all effects connected to
#		# edges that become inactive some time are handled
#		msg.new[which(msg.new<0)] <- msg.new[which(msg.new<0)]+1
#		
#		# which entries became zero? these are the ones becoming active now
#		nowactive <- which(msg>0 & msg.new==0)
#		nowpassive <- which(msg<0 & msg.new==0)
#		coords.nowactive <- cbind((nowactive%%nrow(phi)),(floor(nowactive/nrow(phi)))+1)
#		coords.nowpassive <- cbind((nowpassive%%nrow(phi)),(floor(nowpassive/nrow(phi)))+1)
#		
#		gammasave <- gamma
#		gamma[which(gamma[,tp-1]==1 & gamma[,tp]==0),tp] <- 1
#		
#		## set gamma for endpoints of activated edges
#		if(length(nowactive)>0) {
#			for(i in 1:length(nowactive)) {
#				if(phi[nowactive[i]]==1) {
#					gamma[coords.nowactive[i,2],tp] <- 1
#				}
#				if(phi[nowactive[i]]==2) {
#					gamma[coords.nowactive[i,2],tp] <- 0
#				}
#			}
#		}
#		if(length(nowpassive)>0) {
#			for(i in 1:length(nowpassive)) {
#				if(phi[nowpassive[i]]==1) {
#					gamma[coords.nowpassive[i,2],tp] <- 0
#				}
#			}
#		}
#		
#		## now update the msg.new array
#		# lost activations
#		#        gamma==1 before   gamma=0 now      and activating edge
#		tmp <- which((gamma[,tp-1]==1 & gamma[,tp]==0) & phi==1)
#		targets.coords <- cbind((tmp%%nrow(phi)),(floor(tmp/nrow(phi)))+1)
#		## which ones must be set
#		msg.new[targets.coords[which(msg.new[targets.coords]==0),]] <- -2 # after two timepoints, status will be 0		
#		
#		# activations
#		tmp <- which((phi*gamma[,tp])==1)
#		targets.coords <- cbind((tmp%%nrow(phi)),(floor(tmp/nrow(phi)))+1)
#		## which ones must be set
#		msg.new[targets.coords[which(msg.new[targets.coords]==0),]] <- 1
#		
#		# inhibitions
#		tmp <- which((phi*gamma[,tp])==2)
#		targets.coords <- cbind((tmp%%nrow(phi)),(floor(tmp/nrow(phi)))+1)
#		## which ones must be set
#		msg.new[targets.coords[which(msg.new[targets.coords]==0),]] <- 2		
#		
#		
#		# lost inhibitions
#		# nothing to do here, if inhib active -> status 0, if lost then, status stays
#		# 0 except for another activation comes in, which is handled before
#		# in case of activating edges coming in from active node, reactivate the 
#		# target
#		tmp <- which((phi*gamma[,tp-1])==2 & (phi*gamma[,tp])==0)
#		targets.coords <- cbind((tmp%%nrow(phi)),(floor(tmp/nrow(phi)))+1)
#		## incoming in node to which the ending inhibition points: phi[,targets.coords[2]]
#		## which of these are active now: gamma[,tp]
#		## if one of these active nodes is an activation, set the gamma value to 1 to restore 
#		## the activity
#		if(any(which(phi[,targets.coords[2]] * gamma[,tp]==1))) {
#			gamma[targets.coords[2],tp] <- 1
#		}		
#		tp <- tp + 1
#	}
#	gamma <- gamma[-1,-(maxtp+1)]
#	if(reps>1) {
#		gamma <- expand.gamma(gamma,reps)
#	}
#	return(gamma)
#}

