# TODO: Add comment
# 
# Author: benderc
###############################################################################

## use inhibMCMC output and get mcmc diagnostic, from gelman 2003
gelman_diag <- function(ret) {
	lst <- list()
	for(i in 1:ncol(ret$ltraces)) {
		lst[[i]] <- mcmc(data=ret$ltraces[,i])
	}
	obj <- as.mcmc.list(lst)
	gelman.diag(obj)
}
