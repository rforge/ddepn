# Calculate the Likelihood given an optimal state sequence
# returns the likelihood and parameter point estimates theta
# 
# Author: benderc
###############################################################################
# needs 	library(genefilter)
likl <- function(dat, gammax, scale_lik=FALSE) {
	ieg0 <- ifelse(gammax==0,1,NA)
	ieg1 <- ifelse(gammax==1,1,NA)
	rMeans0=rowMeans(dat*ieg0, na.rm=TRUE)
	rMeans1=rowMeans(dat*ieg1, na.rm=TRUE)
	rSds0=rowSds(dat*ieg0, na.rm=TRUE)
	rSds1=rowSds(dat*ieg1, na.rm=TRUE)
	# penalize the likelihood, small differences in means 
	# give smaller likelihood
	penalty <- -log(abs(rMeans0-rMeans1))
	penalty[is.na(penalty)] <- 0
	#if(shrink) {
	#	L=log(pmin(dnorm(rep(rMeans0,ncol(dat)), mean=rMeans0, sd=rSds0)*ieg0,
	#					dnorm(rep(rMeans1,ncol(dat)), mean=rMeans1, sd=rSds1)*ieg1,na.rm=TRUE)) - penalty
	#} else {
	L=log(pmin(dnorm(dat, mean=rMeans0, sd=rSds0)*ieg0,
					dnorm(dat, mean=rMeans1, sd=rSds1)*ieg1,na.rm=TRUE)) - penalty
	#}
	## perform scaling of the likelihoods according to the number of datapoints that are in each class
	if(scale_lik) {
		scaleones <- matrix(rep(rowSums(gammax), ncol(gammax)),nrow=nrow(gammax), ncol=ncol(gammax))*ieg1
		scalezeros <- matrix(rep(rowSums(1-gammax), ncol(gammax)),nrow=nrow(gammax), ncol=ncol(gammax))*ieg0
		L <- pmin(L/scaleones,L/scalezeros,na.rm=TRUE)
	}
	return(list(L=L, theta=cbind(mu.active=rMeans1, sd.active=rSds1, mu.passive=rMeans0, sd.passive=rSds0)))
}
