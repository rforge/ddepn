# Calling function to invoke the modelling framework
# Data matrix dat in the following format:
# 	rows: antibodies/nodes of network
#	cols: experiments, labeled by EXPLABEL_time
# lambda: prior hyperparameter, try Wehrli/Husmeier sparsity prior, lambda is the same as beta in this
#         paper
# B     : prior weights matrix
# fan.in: maximum number of incoming edges, used for efficient computation of the 
#         normalization factor for the prior distribution
# Author: benderc
###############################################################################

ddepn <- function(dat, phiorig=NULL, phi=NULL, stimuli=NULL, th=0.5, inference="netga", pdf=NULL,
                  multicores=FALSE, maxiterations=1000, p=500, q=0.3, m=0.8, P=NULL,
				  usebics=TRUE, cores=2, lambda=NULL, B=NULL, maxiter=100, fanin=4) {
	# get the experiments, i.e. the stimuli/inhibitor combinations, if not provided
	# works if format of dat is like:
	# colnames contain the experiments in form STIMULUS_time
	if(is.null(stimuli)) {
		cols <- colnames(dat)
		tmp <- sapply(cols, function(x) strsplit(x,"_")[[1]])
		exps <- unique(tmp[1,])
		stims <- sapply(exps, function(x) strsplit(x,"&")[[1]])
		allstim <- unique(unlist(stims))
		stimuli <- list()
		for(i in 1:length(stims)) {
			el <- stims[[i]]
			x <- match(el,allstim)
			names(x) <- el
			stimuli[[i]] <- x
		}
	} else {
		allstim <- unique(names(unlist(stimuli)))	
	}
	# test if the stimuli are included as nodes, if not, add them
	test <- which(is.na(match(allstim,rownames(dat))))
	if(length(test)>0) {
		for(tt in test) {
			dat <- rbind(rep(0.0,ncol(dat)), dat)
			rownames(dat) <- c(allstim[tt],rownames(dat))
		}
	}
	n <- nrow(dat)
	phinames <- rownames(dat)
	
	## if information for prior is given then test if it is complete and
	## compute the normalisation factor Z(lambda) as in Fröhlich2007 or Wehrli/Husmeier 2007
	if(!is.null(lambda) || !is.null(B)) {
		if(is.null(lambda) || is.null(B)) {
			print("Prior information incomplete. Please provide arguments lambda and B.")
			stop()
		}
		## get normalisation factor for the networks
		print("Computing prior normalisation factor...")
		Z <- zlambda(B, lambda)
		print("done.")
	} else {
		Z <- NULL
	}
	
	if(inference=="netga") {
		scorefile <- paste("score",sub("\\.pdf","",pdf),".pdf",sep="")
        stime <- system.time(P <- netga(dat,stimuli,P=P,maxiterations=maxiterations,p=p,q=q,m=m,multicores=multicores,usebics=usebics,cores=cores,lambda=lambda,B=B,Z=Z,maxiter=maxiter,scorefile=scorefile,fanin=fanin))
        #if(is.null(phiorig)) {
        #	phiorig <- matrix(0,nrow=nrow(dat),ncol=nrow(dat),dimnames=list(rownames(dat),rownames(dat)))
        #}
        phi.activation.count <- phi.inhibition.count <- weights.tc <- matrix(0,nrow=nrow(dat),ncol=nrow(dat),dimnames=list(rownames(dat),rownames(dat))) 
		# now compare all graphs to the original
        result <- NULL
		resultep <- NULL
        for(i in 1:length(P)) {
			if(!is.null(phiorig))
            	result <- rbind(result, cbind(compare.graphs.tc(phiorig, P[[i]]$phi),t(as.matrix(stime))))
			#resultep <- rbind(resultep, cbind(compare.graphs.ep(phiorig, P[[i]]$phi,tps=1:50,stimuli=stimuli), t(as.matrix(stime))))
			phi.activation.count[which(P[[i]]$phi==1)] <- phi.activation.count[which(P[[i]]$phi==1)] + 1
			phi.inhibition.count[which(P[[i]]$phi==2)] <- phi.inhibition.count[which(P[[i]]$phi==2)] + 1
			weights.tc <- weights.tc + detailed.to.simple.regulations(P[[i]]$phi)
		}
        conf.act <- round(phi.activation.count/length(P),digits=3)
        conf.inh <- round(phi.inhibition.count/length(P),digits=3)
        weights.tc <- round(weights.tc/length(P),digits=3)
        ret <- list(dat=dat, phi.activation.count=phi.activation.count,
					phi.inhibition.count=phi.inhibition.count,
					phi.orig=phiorig, phi=NULL, weights=NULL,
					weights.tc=weights.tc, result=result, conf.act=conf.act,conf.inh=conf.inh,
					stimuli=stimuli) #, resultep=resultep)
        ret <- get.phi.final(ret,th)
        plotrepresult(ret,pdf)
        ret[["P"]] <- P
	} else {
		if(inference=="mcmc") {
			#phistart <- matrix(sample(c(0,1,2),n*n,replace=T), nrow=n, ncol=n, dimnames=list(phinames,phinames))
			phistart <- matrix(0, nrow=n, ncol=n, dimnames=list(phinames,phinames))
			ret <- mcmc_ddepn(dat, phiorig=phiorig, phi=phistart, stimuli=stimuli,
						th=th, multicores=multicores, pdf=pdf, maxiterations=maxiterations,
						usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, maxiter=maxiter,fanin=fanin)
		}
	}
	ret
}

