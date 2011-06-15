# TODO: Add comment
# 
# Author: benderc
###############################################################################
calcpr <- function(lambda, B, phi, gam) {
	-log(2) - log(lambda) + (-(abs(B - phi)^gam))/lambda
}

prior <- function(phi, lambda=NULL, B=NULL, Z=NULL, gam=NULL, it=NULL, K=NULL, priortype="laplaceinhib") {
	if(priortype=="laplace") {
		## if B only holds the information if there is an edge or not, then use laplace. if 
		## additionally the type of edge is encoded in B, then use laplaceinhib
		phi <- detailed.to.simple.regulations(phi)
		if(is.null(gam))
			gam <- 1
		## integrate if lambda==NA
		## integration should be over lambda in [0,Inf], since I am integrating over the logged 
		## version, use lambda in [1, Inf]
		## don't know if this is correct, integral to Inf is not bounded, and doesn't work
		## so I define a certain range for lambda (1:100)
		if(is.na(lambda)) {
			PGlambda <- sum(sapply(1:length(B), function(i,B,phi,gam) integrate(calcpr, 1,100, B=B[i], phi=phi[i], gam=gam)$value, B=B, phi=phi, gam=gam))			
		} else {
			PGlambda <- sum(calcpr(lambda, B, phi, gam))
		}
		#PGlambda <- sum(-log(2) - log(lambda) + (-(abs(B - phi)^gam))/lambda)
	} else if(priortype=="scalefree") {
		PGlambda <- pgs(phi,gam,K,it) ## logged value is returned
	} else if(priortype=="laplaceinhib") {
		# B has same dimensions as phi, reduce phi to network with only one edge type
		## make sure that the difference of an inhibition edge in stead of no edge
		## is the same as an activation edge in stead of no edge
		## the difference of an activation in stead of an inhibition or vice
		## versa should be high, since introducing the wrong effect is worse than
		## leaving out the effect
		if(is.null(gam))
			gam <- 1
		phi[phi==2] <- -1	
		## integrate if lambda==NA
		## integration should be over lambda in [0,Inf], since I am integrating over the logged 
		## version, use lambda in [1, Inf]
		## don't know if this is correct, integral to Inf is not bounded, and doesn't work
		## so I define a certain range for lambda
		if(is.na(lambda)) {
			PGlambda <- sum(sapply(1:length(B), function(i,B,phi,gam) integrate(calcpr, 1,100, B=B[i], phi=phi[i], gam=gam)$value, B=B, phi=phi, gam=gam))
		} else {
			PGlambda <- sum(calcpr(lambda, B, phi, gam))
		}
		#PGlambda <- sum(-log(2) - log(lambda) + (-(abs(B - phi)^gam))/lambda)
	} else {
		## uniform prior: set prior prob to 1 for all possible networks,
		## will add a log-prior of 0 to the likelihood to obtain the 
		## posterior. Thus, equals simple likelihood optimisation
		PGlambda <- log(1)
	}
	PGlambda
}
