# Calculate the Likelihood given an optimal state sequence
# returns the likelihood and parameter point estimates theta
# 
# Author: benderc
###############################################################################
# needs 	library(genefilter)
likl <- function(dat, gammax) {
	rMeans0=rowMeans(dat*ifelse(gammax==0,1,NA), na.rm=TRUE)
	rMeans1=rowMeans(dat*ifelse(gammax==1,1,NA), na.rm=TRUE)
	rSds0=rowSds(dat*ifelse(gammax==0,1,NA), na.rm=TRUE)
	rSds1=rowSds(dat*ifelse(gammax==1,1,NA), na.rm=TRUE)
	L=log(pmin(dnorm(dat, mean=rMeans0, sd=rSds0)*ifelse(gammax==0,1,NA),
				dnorm(dat, mean=rMeans1, sd=rSds1)*ifelse(gammax==1,1,NA),na.rm=TRUE))
	## perform scaling of the likelihoods according to the number of datapoints that are in each class
	scaleones <- matrix(rep(rowSums(gammax), ncol(gammax)),nrow=nrow(gammax), ncol=ncol(gammax))*ifelse(gammax==1,1,NA)
	scalezeros <- matrix(rep(rowSums(1-gammax), ncol(gammax)),nrow=nrow(gammax), ncol=ncol(gammax))*ifelse(gammax==0,1,NA)
	L <- pmin(L/scaleones,L/scalezeros,na.rm=TRUE)
  	return(list(L=L, theta=cbind(mu.active=rMeans1, sd.active=rSds1, mu.passive=rMeans0, sd.passive=rSds0)))
}
