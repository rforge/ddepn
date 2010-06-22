# Calculate the Likelihood of datx given gammax and thetax
# 
# Author: benderc
###############################################################################
# needs 	library(genefilter)
likl <- function(dat, gammax) {
	rMeans0=rowMeans(dat*ifelse(gammax==0,1,NA), na.rm=TRUE)
	rMeans1=rowMeans(dat*ifelse(gammax==1,1,NA), na.rm=TRUE)
	rSds0=rowSds(dat*ifelse(gammax==0,1,NA), na.rm=TRUE)
	rSds1=rowSds(dat*ifelse(gammax==1,1,NA), na.rm=TRUE)
	L=log2(pmin(dnorm(dat, mean=rMeans0, sd=rSds0)*ifelse(gammax==0,1,NA),
				dnorm(dat, mean=rMeans1, sd=rSds1)*ifelse(gammax==1,1,NA),na.rm=TRUE))
  	return(list(L=L, theta=cbind(mu.active=rMeans1, sd.active=rSds1, mu.passive=rMeans0, sd.passive=rSds0)))
}
