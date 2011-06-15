# For inhibMCMC results: use multiple MCMC chains and retrieve significant edges
# 
# Author: benderc
###############################################################################

create_signetwork <- function(ret, alpha=0.001, adj.method="BY", plot=FALSE, type="wilcox", alternative="one.sided", paired=FALSE, ord=NULL) {
	## check if only the samplings list is given, or the whole mcmc return object
	if(!is.null(ret$samplings)) {
		ret <- ret$samplings
	}
	if(length(ret)<=2) {
		stop(paste("Caution: Only", length(ret), "parallel MCMC chain(s) found and statistical testing for", 
				   "edges not possible. Please run at least 3 parallel chains."))
	}
	prior <- ret[[1]]$phi.orig
	stimuli <- ret[[1]]$stimuli
	## use statistical testing to get the consensus network
	if(is.null(ord))
		ord <- rownames(ret[[1]]$phi)
	gsnt <- perform_network_tests(ret, alpha=alpha, plot=plot, type=type, alternative=alternative, paired=paired, ord=ord)
	pacts <- p.adjust(gsnt$netact,method=adj.method)
	pinhs <- p.adjust(gsnt$netinh,method=adj.method)
	#inferred <- matrix(0, nrow=nrow(pacts), ncol=ncol(pacts), dimnames=dimnames(pacts))
	inferred <- matrix(0, nrow=nrow(gsnt$network), ncol=ncol(gsnt$network), dimnames=dimnames(gsnt$network))
	acts <- which(pacts<=alpha)
	inhs <- which(pinhs<=alpha)
	inferred[setdiff(acts, inhs)] <- 1
	inferred[setdiff(inhs, acts)] <- 2
	inferred
}
## calls perform_edge_test for each pair of from/to nodes
#
#ord=c("EGF","H","P","E","HRG","AKT","ERBB1","ERBB2","ERBB3","ERK12","GSK3","MEK12", 
#		"MTOR","p38","p70S6K","PDK1","PKC","PLCG","PRAS","cSRC"),
perform_network_tests <- function(ret, alpha=0.05, plot=FALSE, type="wilcox", alternative="one.sided", paired=FALSE, ord=NULL) {
	## get protein names and order, if wanted 
	pA <- rownames(ret[[1]]$phi.orig)
	if(!is.null(ord)) {
		dnm <- match(ord,pA)
		dnm <- dnm[!is.na(dnm)]
		dnames <- pA[dnm]
	} else {
		dnames <- pA
	}
	## define adjacency matrix
	network <- netact <- netinh <- matrix(0, nrow=length(dnames), ncol=length(dnames), dimnames=list(dnames,dnames))
	stimuli <- unique(names(c(unlist(ret[[1]]$stimuli))))
	## check each edge
	for(ai in 1:length(dnames)) {
		print(dnames[ai])
		for(bi in 1:length(dnames)) {
			if(ai==bi || bi %in% stimuli)
				next
			from <- dnames[ai]
			to <- dnames[bi]
			rr <- perform_edge_test(ret, from, to, alpha=alpha, plot=plot, type=type, alternative=alternative, paired=paired)
			netact[from,to] <- rr[1]
			netinh[from,to] <- rr[2]
			if(rr[1]<=alpha) {
				network[from,to] <- 1
			}
			if(rr[2]<=alpha) {
				network[from,to] <- 2
			}
		}
	}
	list(network=network, netact=netact, netinh=netinh)
}
## calculate tests if activation of from is significantly higher than activation of to
## same for inhibitions
perform_edge_test <- function(ret, from, to, alpha=0.05, plot=TRUE, type="wilcox", alternative="one.sided", paired=FALSE) {
	## get activation, inhibition counts
	aact <- sapply(ret, function(x) x$conf.act[from,to])
	ainh <- sapply(ret, function(x) x$conf.inh[from,to])
	if(alternative=="two.sided") {
		if(type=="ttest") {
			ag <- signif(t.test(aact,ainh,alternative="two.sided",paired=paired)$p.value,digits=3)
			al <- median(aact)>median(ainh)
		} else {
			ag <- signif(wilcox.test(aact,ainh,alternative="two.sided",paired=paired)$p.value,digits=3)
			al <- median(aact)>median(ainh)
		}	
	} else {
		if(type=="ttest") {
			ag <- signif(t.test(aact,ainh,alternative="greater",paired=paired)$p.value,digits=3)
			al <- signif(t.test(aact,ainh,alternative="less",paired=paired)$p.value,digits=3)
		} else {
			ag <- signif(wilcox.test(aact,ainh,alternative="greater",paired=paired)$p.value,digits=3)
			al <- signif(wilcox.test(aact,ainh,alternative="less",paired=paired)$p.value,digits=3)
		}	
	}
	mat <- cbind(aact,ainh)
	colnames(mat) <- c("activation","inhibition")
	main <- c(paste(from,"-",to), type)#,
		#paste("p(act>inh) =",ag),
		#paste("p(act<inh) =",al))
	if(plot) {
		boxplot(as.data.frame(mat),main=main,
				las=1, ylim=c(0,1.2), ylab="edge confidence")
		tx <- c(a_greater=ag,a_less=al)
		draw_segments(tx, range(mat), alpha)
	}
	invisible(c(ag=ag, al=al))
}
## draw significance bars over boxplots
# comes from here:
#  source("~/projects/rfunctions/datavisualisation/draw_segments.R")
draw_segments <- function(tx, rng, alpha) {
	## which pvalues are significant, alternative greater or less?
	wt <- which(tx<=alpha)
	lng <- 0.1
	if(length(wt)>0) {
		#browser()
		for(i in 1:length(wt)) {
			wti <- wt[i]
			txi <- tx[wti]
			wti <- ifelse(wti<3,1,2)
			#wti <- wti %% 2 +1
			x0 <- wti*2-1
			x1 <- wti*2
			y0 <- rng[2] + 2*lng
			y1 <- rng[2] + 2*lng
			segments(x0,y0,x1,y1)
			y0 <- rng[2] + 2*lng#-lng*rng[2]
			y1 <- rng[2] + lng
			segments(x0,y0,x0,y1)
			segments(x1,y0,x1,y1)
			sigt <- "*"
			if(txi <= 0.01)
				sigt <- "**"
			if(txi <= 0.001)
				sigt <- "***"
			altern <- switch(names(wti),
					a_greater=altern<-"(act > inh)",
					a_less=altern<-"(act < inh)",
					b_greater=altern<-"(act > inh)",
					b_less=altern<-"(act < inh)")
			#txt <- paste(sigt, " (", expression(paste(p[greater])),"=",txi,")",sep="")
			#txt <- expression(paste(sigt, " (p[greater]=", txi, ")",sep=""))
			txt <- paste(sigt, " (p=",txi,")",sep="")
			text(x0+0.5, (rng[2]+2.5*lng), txt)
			text(x0+0.5, (rng[2]+1.5*lng), altern)
		}
	}
}
