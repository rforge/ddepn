# Produce ROC curve and calculate AUC for MCMC intermediate result
# 
# Author: benderc
###############################################################################
trapezoid <- function(sp,sn) {
	sum(diff(sp)*(sn[-1]+sn[-length(sn)]))/2
}

mcmc_performance <- function(lst, plot=TRUE) {
	thlim <- seq(0,1,by=0.01)
	phiorig <- lst$phi.orig
	stats <- NULL
	for(th in thlim) {
		lst <- get.phi.final(lst,th=th)
		if(is.null(phiorig)) {
			comp <- rep(0,8)
			names(comp) <- c("tp","tn","fp","fn","sn","sp","prec","f1")
		} else {
			comp <- compare.graphs.tc(phiorig=phiorig,phi=lst$phi,ignore.type=FALSE)
		}
		stats <- rbind(stats,comp[1:6])
	}
	x <- c(0,stats[,"sp"],1)
	y <- c(1,stats[,"sn"],0)
	auc <- trapezoid(x,y)
	if(plot) {
		plot((1-x), y, type="l", lty=1, #"dashed",
			xlim=c(0,1),ylim=c(0,1), col="black", xlab="1-SP", ylab="SN",
			main=c("ROC curve",paste("AUC: ",signif(auc,digits=2),sep="")),lwd=2)
	}
	list(stats=stats,auc=auc)
}


