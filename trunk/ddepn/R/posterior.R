# Calculate posterior probability using prior as in Fröhlich 2007/Wehrli/Husmeier 2007
# or by sparsity prior
# 
# Author: benderc
###############################################################################

# scalefree prior
#posterior <- function(phi, L, lambda, B=NULL, Z=NULL, gam=2.2, it=500, K=0.8) {
posterior <- function(phi, L, lambda=NULL, B=NULL, Z=NULL, gam=NULL, it=NULL, K=NULL, priortype="laplaceinhib") {
	if(priortype %in% c("laplaceinhib","laplace","scalefree")) {
		PGlambda <- prior(phi, lambda, B, Z, gam, it, K, priortype)
		post <- L + PGlambda
	#laplace <- !is.null(lambda) && !is.null(B) && !is.null(Z)
	#scalefree <- !is.null(gam) && !is.null(it) && !is.null(K)
#	if(priortype=="laplace") {
#		EG <- -(abs(B - detailed.to.simple.regulations(phi))/lambda)
#		prefix <- -log(2) - log(lambda)
#		PGlambda <- sum(prefix + EG)
#		####EG <- sum(abs(B - detailed.to.simple.regulations(phi)))
#		#####PGlambda <- log2(exp(-lambda * EG)) - Z
#		#####PGlambda <- log2(2^(-lambda * EG)) - Z
#		####PGlambda <- (-lambda * EG) - Z
#		post <- L + PGlambda
#	} else if(priortype=="scalefree") {
#		post <- L + log(pgs(phi,gam,K,it))
#	} else if(priortype=="laplaceinhib") {
#		# B has same dimensions as phi, reduce phi to network with only one edge type
#		## make sure that the difference of an inhibition edge in stead of no edge
#		## is the same as an activation edge in stead of no edge
#		## the difference of an activation in stead of an inhibition or vice
#		## versa should be high, since introducing the wrong effect is worse than
#		## leaving out the effect
#		phi[phi==2] <- -1	
#		PGlambda <- sum(-log(2) + log(lambda) + (-abs(B - phi))/lambda)
#		post <- L + PGlambda
	} else {
		stop("posterior.R: Error - Prior information not specified correctly.")
	}
	post
}


pi <- function(i,gam,N) {
	(1-(1/(1-gam)))/(N^(1-(1/(1-gam)))) * i^(-(1/(gam-1)))
}

## phi: network
## N : number of nodes
## gam: degree distribution coefficient: P(K) ~ K^gam
## perm: permutation
pgs <- function(phi,gam,K=0.8,it=500) {
	N <- nrow(phi)
	ind <- which(phi!=0,arr.ind=TRUE)
	fixe <- exp(-N*K*(1-sum(pi(1:N,gam,N)^2)))
	res <- 0
	piN <- pi(1:N, gam, N)
	mu <- 1/(gam-1)
	for(b in 1:it) {
		perm <- sample(1:N)
		piNp <- piN[perm]
		res <- res + prod(exp(2*N*K*piNp[ind[,1]]*piNp[ind[,2]]-1)) * fixe
	}
	res <- res/it
	res
}
