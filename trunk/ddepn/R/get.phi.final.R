# TODO: Add comment
# 
# Author: benderc
###############################################################################


# get the final network
get.phi.final <- function(lst,th=0.8) {
	n <- nrow(lst$dat)
	rn <- rownames(lst$dat)
	phi.final <- matrix(0,nrow=n,ncol=n,dimnames=list(rn,rn))
	weights <- matrix(0,nrow=n,ncol=n,dimnames=list(rn,rn))
	#weights.tc <- lst$weights.tc
	#phi.activation.count <- lst$phi.activation.count
	#phi.inhibition.count <- lst$phi.inhibition.count
	ca <- signif(lst$conf.act,digits=3)
	ci <- signif(lst$conf.inh,digits=3)
	for(i in 1:nrow(phi.final)) {
		for(j in 1:ncol(phi.final)) {
			if(i==j) {
				weights[i,j] <- 0
				phi.final[i,j] <- 0
			} else {
				if(ca[i,j]>=th || ci[i,j]>=th) {
					if(ca[i,j]==ci[i,j]) {
						phi.final[i,j] <- 1
						#weights[i,j] <- paste(ca[i,j],"??")
						weights[i,j] <- "undisting"
					}
					if(ca[i,j]>ci[i,j]) {
						phi.final[i,j] <- 1
						weights[i,j] <- ca[i,j]
					}
					if(ca[i,j]<ci[i,j]) {
						phi.final[i,j] <- 2
						weights[i,j] <- ci[i,j]
					}
				}
			}
		}
	}
	lst$phi <- phi.final
	lst$weights <- weights
	return(lst)
}

