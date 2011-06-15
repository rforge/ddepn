# TODO: Add comment
# 
# Author: benderc
###############################################################################
get.phi.final.mcmc <- function(retlist,maxiterations,prob=.333,qu=.99999) {
	## find significant edges
	## get .9999 quantile for the number of successes x with P(x<=X) >= .9999, d.h. P(x>X)<.0001
	## for p = 0.33. This corresponds to random drawing of activating/inhibiting edges. If a higher
	## number of edges is found in the network, include the edge
	#prob <- .3333
	cutrng <- c(qbinom(qu,maxiterations,prob,lower.tail=TRUE),
			qbinom(qu,maxiterations,prob,lower.tail=FALSE))
	for(rtl in 1:length(retlist)) {
		retX <- retlist[[rtl]]
		acts <- which(retX$freqa>=cutrng[2])
		inhs <- which(retX$freqi>=cutrng[2])
		aimax <- apply(cbind(retX$freqa[intersect(acts,inhs)],retX$freqi[intersect(acts,inhs)]),1,function(x) which(x==max(x))[1]) ## when the same number of inhs and acts, take act
		intai <- intersect(acts,inhs)
		acts <- acts[-match(intai[which(aimax==2)],acts)]
		inhs <- inhs[-match(intai[which(aimax==1)],inhs)]
		net <- retX$phi
		net[net!=0] <- 0
		weights <- net
		net[acts] <- 1
		weights[acts] <- retX$conf.a[acts]
		net[inhs] <- 2
		weights[inhs] <- retX$conf.i[inhs]
		weights <- signif(weights,digits=3)
		retX$phi <- net
		retX$weights <- weights
		retlist[[rtl]] <- retX
		#plotdetailed(net)
	}
	retlist
}

# get the final network
get.phi.final <- function(lst,th=0.8) {
	n <- nrow(lst$dat)
	rn <- rownames(lst$dat)
	phi.final <- matrix(0,nrow=n,ncol=n,dimnames=list(rn,rn))
	weights <- matrix(0,nrow=n,ncol=n,dimnames=list(rn,rn))
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

