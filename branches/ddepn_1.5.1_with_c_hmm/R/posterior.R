# Calculate posterior probability using prior as in Fröhlich 2007/Wehrli/Husmeier 2007
# or by sparsity prior
# 
# Author: benderc
###############################################################################

# scalefree prior
#posterior <- function(phi, L, lambda, B=NULL, Z=NULL, gam=2.2, it=500, K=0.8) {
posterior <- function(phi, L, lambda=NULL, B=NULL, Z=NULL, gam=NULL, it=NULL, K=NULL, priortype="laplaceinhib") {
	if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {
		PGlambda <- prior(phi, lambda, B, Z, gam, it, K, priortype)
		post <- L + PGlambda
	} else {
		stop("posterior.R: Error - Prior information not specified correctly.")
	}
	post
}

## exact formula
pri <- function(i,gam=2.2,N) {
	mu <- 1/(gam-1)
	den <- 0
	for(j in 1:N) {
		den <- den + j^-mu
	}
	prival <- (i^-mu) / den
	prival
}

## phi: network
## N : number of nodes
## gam: degree distribution coefficient: P(K) ~ K^gam
## perm: permutation
pgs <- function(phi,gam,K=0.8,it=500) {
	N <- nrow(phi)
	ind <- which(phi!=0,arr.ind=TRUE) ## edge indices
	ind2 <- which(phi==0,arr.ind=TRUE) ## indices of not connected node pairs
	piN <- pri(1:N, gam, N) ## probabilites of nodes i
	res <- 0 ## sum of probs over iterations
	for(b in 1:it) {
		perm <- sample(1:N) ## permutation
		piNp <- piN[perm] ## probability permutation
		e1 <- log(1 - exp(-2*N*K*piNp[ind[,1]]*piNp[ind[,2]])) ## edge probabilities
		e2 <- -2*N*K*piNp[ind2[,1]]*piNp[ind2[,2]] ## missing edge probabilities
		res <- res + exp(sum(e1)+sum(e2))
	}
	res <- log(res) - log(it)
	res
} 
