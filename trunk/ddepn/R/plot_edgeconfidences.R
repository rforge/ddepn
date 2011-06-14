# For an inhibMCMC result object, plot the confidences of 
# activations and inhibitions for each protein
# 
# Author: benderc
###############################################################################


make_edge_df <- function(samp) {
	sam <- samp[[1]]	
	nr <- nrow(sam)
	L <- length(samp)
	## sequence of source nodes for all 10 runs
	src <- rep(rep(rownames(sam),nr),L)
	src <- factor(src,levels=rownames(sam))
	## sequence of dest nodes for all 10 runs
	dest <- rep(rep(rownames(sam),each=nr),L)
	dest <- factor(dest,levels=rownames(sam))
	## number of run
	nr <- rep(1:L,each=nr*nr)
	## counts
	y <- sapply(samp, "[")
	df <- data.frame(y=as.vector(y), src=src, dest=dest)
	df
}

plot_edgeconfidences <- function(ret, start=1, stop=NULL, act="conf.act", inh="conf.inh",...) {
	#stopifnot(require(lattice))
	## check if only the samplings list is given, or the whole mcmc return object
	if(!is.null(ret$samplings)) {
		ret <- ret$samplings
	}
	## get the confidences/counts for each edge
	sampa <- lapply(ret, function(x) x[[act]])
	sampi <- lapply(ret, function(x) x[[inh]])
	## transform to data frame for lattice plotting
	dfa <- make_edge_df(sampa)
	dfi <- make_edge_df(sampi)
	type <- paste(c(as.character(dfa[,"src"]),as.character(dfi[,"src"])),c(rep("+",nrow(dfa)),rep("-",nrow(dfi))),sep="_")
	df <- cbind(rbind(dfa, dfi),type)
	legendX <- simpleKey(c("activation (+)","inhibition (-)"),points=FALSE,rectangles=TRUE,lines=FALSE,coL=c("#0000FF","#FF0000"))
	legendX$corner <- c(1,1)
	legendX$x <- 0.9
	legendX$y <- 0.9
	legendX$rectangles$col <- c("#0000FF","#FF0000")
	legendX$text$cex <- 2
	## make a lattice boxplot figure
	bwplot(y ~ type | dest, data=df, pch=NULL, main=list(label="Confidences of activation/inhibition edges over 10 MCMC runs (counted after burnin phase)",cex=2),
			ylab=list(label="confidence",cex=2),xlab=list(label="source node", cex=2),
			scales=list(rot=90,alternating=3,cex=1.5),
			panel=function(...) {
				panel.grid(h = 0, v = 17, col="#111111aa", lty=3, lwd=2)
				panel.bwplot(...,fill=rep(c("blue","red"),18))
			},
			key=legendX,
			par.settings=list(layout.heights=list(bottom.padding=10)))
}

