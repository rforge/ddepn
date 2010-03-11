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
