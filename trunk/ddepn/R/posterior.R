# Calculate posterior probability using prior as in Wehrli/Husmeier 2007
# 
# Author: benderc
###############################################################################


posterior <- function(phi, L, lambda, B, Z) {
	# B has same dimensions as phi, reduce phi to network with only one edge type
	EG <- sum(abs(B - detailed.to.simple.regulations(phi)))
	#PGlambda <- log2(exp(-lambda * EG)) - Z
	#PGlambda <- log2(2^(-lambda * EG)) - Z
	PGlambda <- (-lambda * EG) - Z
	post <- L + PGlambda
	post
}


posterior2 <- function(phi, L, lambda, B=NULL, Z=NULL, gam=2.2, it=500, K=0.8) {
	post <- L + log2(pgs(phi,gam,K,it))
	post
}

pi <- function(i,gam,N) {
	(1-(1/(1-gam)))/(N^(1-(1/(1-gam)))) * i^(-(1/(gam-1)))
}

## phi: network
## N : number of nodes
## gam: degree distribution coefficient: P(k) ~ k^gam
## perm: permutation
pgs <- function(phi,gam,K=0.8,B=500) {
	N <- nrow(phi)
	ind <- which(phi!=0,arr.ind=TRUE)
	fixe <- exp(-N*K*(1-sum(pi(1:N,gam,N)^2)))
	res <- 0
	piN <- pi(1:N, gam, N)
	mu <- 1/(gam-1)
	for(b in 1:B) {
		perm <- sample(1:N)
		piNp <- piN[perm]
		res <- res + prod(exp(2*N*K*piNp[ind[,1]]*piNp[ind[,2]]-1)) * fixe
	}
#	for(b in 1:B) {
#		perm <- cbind(sample(1:5),1:5)
#		ind2 <- ind
#		ind2[,1] <- perm[match(ind[,1],perm[,2]),1]
#		ind2[,2] <- perm[match(ind[,2],perm[,2]),1]
#
#		piind <- pi(ind2,gam,N)
#		res <- res + prod(exp(2*N*K*piind[,1]*piind[,2]-1)) * fixe
#	}
	res <- res/B
	res
}
