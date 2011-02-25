# Compute normalisation factor for prior as defined in Wehrli/Husmeier 2007.

# Z = prod_v ( sum_pa(v) ( exp(-lambda * (sum_minpa(v)(1-B[v,m]) + sum_mnotinpa(v)(B[v,m])) ) ) )
#   =                                           B[-v,pa(v)]
# Author: benderc
###############################################################################

zlambda <- function(B, lambda) {
	V <- 1:nrow(B) #rownames(B)
	Z <- 0
	for(v in 1:nrow(B)) {
		print(v)
		# get parental configurations for maximum 3 incoming edges (fan.in always 3)
		## one parent
		Vp <- V[-v]
		summe <- 0
		for(i in 1:length(Vp)) {
			# parents of v
			pa <- Vp[i]
			summand <- exp(-lambda * (sum(1-B[pa,v]) + sum(B[setdiff(Vp,pa),v])))
			summe <- summe + summand
		}
		## two parents
		for(i in 1:length(Vp)) {
			for(j in (i+1):length(Vp)) {
				if(j>length(Vp))
					next
				pa <- Vp[c(i,j)]
				summand <- exp(-lambda * (sum(1-B[pa,v]) + sum(B[setdiff(Vp,pa),v])))
				summe <- summe + summand
			}
		}
		## three parents
		for(i in 1:length(Vp)) {
			for(j in (i+1):length(Vp)) {
				if(j>length(Vp))
					next
				for(k in (j+1):length(Vp)) {
					if(k>length(Vp))
						next
					pa <- Vp[c(i,j,k)]
					summand <- exp(-lambda * (sum(1-B[pa,v]) + sum(B[setdiff(Vp,pa),v])))
					summe <- summe + summand
				}
			}
		}
		Z <- Z + log2(summe)
	}
	return(Z)
}

