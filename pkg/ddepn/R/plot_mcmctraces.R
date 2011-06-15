# Plot MCMC traces from inhibMCMC chains
# requires package coda
#
# Author: benderc
###############################################################################


plot_mcmctraces <- function(ret, thin=1) {
	samp <- ret$samplings
	ltr <- ret$ltraces
	lst <- list()
	for(i in 1:ncol(ltr)) {
		mcmcobj <- mcmc(ltr[,i],thin=thin)
		lst[[i]] <- mcmcobj
	}
	mcmclst <- mcmc.list(lst)
	plot(mcmclst, main="Posterior traces, 10 chains",auto.layout=FALSE,density=FALSE)
}
