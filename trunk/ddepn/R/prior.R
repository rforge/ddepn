# TODO: Add comment
# 
# Author: benderc
###############################################################################

prior <- function(phi, lambda=NULL, B=NULL, Z=NULL, gam=NULL, it=NULL, K=NULL, priortype="laplaceinhib") {
	if(priortype=="laplace") {
		## if B only holds the information if there is an edge or not, then use laplace. if 
		## additionally the type of edge is encoded in B, then use laplaceinhib
		phi <- detailed.to.simple.regulations(phi)
		PGlambda <- sum(-log(2) - log(lambda) + (-(abs(B - phi)^gam))/lambda)
	} else if(priortype=="scalefree") {
		PGlambda <- log(pgs(phi,gam,K,it))
	} else if(priortype=="laplaceinhib") {
		# B has same dimensions as phi, reduce phi to network with only one edge type
		## make sure that the difference of an inhibition edge in stead of no edge
		## is the same as an activation edge in stead of no edge
		## the difference of an activation in stead of an inhibition or vice
		## versa should be high, since introducing the wrong effect is worse than
		## leaving out the effect
		if(is.null(gam))
			gam <- 2
		phi[phi==2] <- -1	
		PGlambda <- sum(-log(2) - log(lambda) + (-(abs(B - phi)^gam))/lambda)
	}
	PGlambda
}
