# Calling function to invoke the modelling framework
# Data matrix dat in the following format:
# 	rows: antibodies/nodes of network
#	cols: experiments, labeled by EXPLABEL_time
# lambda: laplace prior hyperparameter
# B     : prior weights matrix
# gam   : sparsity prior hyperparameter
# it    : sparsity prior hyperparameter
# k     : sparsity prior hyperparameter
# fan.in: maximum number of incoming edges, used for efficient computation of the 
#         normalization factor for the prior distribution
#
#
#  Priors: laplace: Fröhlich 2007 / Wehrli/Husmeier 2007
#          sparsity: Kamimura 200??? F7 414, Lee 2005 etc., not sure about the correct citation
# Author: benderc
###############################################################################

ddepn <- function(dat, phiorig=NULL, phi=NULL, stimuli=NULL, th=0.5, inference="netga", pdf=NULL,
                  multicores=FALSE, maxiterations=1000, p=500, q=0.3, m=0.8, P=NULL,
				  usebics=TRUE, cores=2, 
				  lambda=NULL, B=NULL, samplelambda=TRUE,
				  maxiter=100, fanin=4,
				  gam=NULL,it=NULL,K=NULL,quantL=.5,quantBIC=.5,
				  debug=FALSE) {
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
	
	laplace <- !is.null(lambda) && !is.null(B)
	sparsity <- !is.null(gam) && !is.null(it) && !is.null(K)
	none <- is.null(lambda) && is.null(B) && is.null(gam) && is.null(it) && is.null(K)
	
	if(laplace) {
		## get normalisation factor for the networks
		print("Computing prior normalisation factor...")
		Z <- zlambda(B, lambda)
		print("done.")
	} else if(sparsity || none) {
		Z <- NULL
	} else {
		stop("Error in function arguments. Please specifiy either lambda/gamma for laplace prior, gam/it/K for sparsity prior or none if no prior distribution should be used.")
	}
	
	if(inference=="netga") {
		scorefile <- paste("score",sub("\\.pdf","",pdf),".pdf",sep="")
        stime <- system.time(P <- netga(dat,stimuli,P=P,maxiterations=maxiterations,
						p=p,q=q,m=m,multicores=multicores,usebics=usebics,
						cores=cores,lambda=lambda,B=B,Z=Z,maxiter=maxiter,
						scorefile=scorefile,fanin=fanin,
						gam=gam,it=it,K=K,quantL=quantL,quantBIC=quantBIC))
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
			#browser()
			if(multicores) {
				# start an mcmc run for each core
				# liste, die die dateinamen und die startnetze enthält
				if(is.null(P)) {
					P <- list()
					phistart <- matrix(0, nrow=n, ncol=n, dimnames=list(phinames,phinames))
					for(cr in 1:cores) {
						P[[cr]] <- list(phi=phistart)
					}
				}
				X <- list()
				for(cr in 1:cores) {
					filename <- paste(sub("\\.pdf","",pdf),"_",cr,".pdf",sep="")
					X[[cr]] <- list(phi=P[[cr]]$phi, pdf=filename)
				}
				retlist <- mclapply(X, runmcmc, dat=dat, phiorig=phiorig, phi=phistart, stimuli=stimuli,
						th=th, multicores=multicores, pdf=pdf, maxiterations=maxiterations,
						usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, samplelambda=samplelambda,
						maxiter=maxiter,fanin=fanin, gam=gam, it=it, K=K,
						mc.preschedule=FALSE,mc.cores=cores)
				#runmcmc(X[[1]],dat=dat, phiorig=phiorig, phi=phistart, stimuli=stimuli,
				#		th=th, multicores=multicores, pdf=pdf, maxiterations=maxiterations,
				#		usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, samplelambda=samplelambda,
				#		maxiter=maxiter,fanin=fanin, gam=gam, it=it, K=K)
				### experimental
				if(debug)
					browser()
				retlist <- get.phi.final.mcmc(retlist, maxiterations, prob=.3333, qu=.99999)
				##### end experimental
				## get the likelihood traces
				ltraces <- sapply(retlist,function(x) x$stats[,"MAP"])
				## make sure ltraces has less than 10000 rows
				if(nrow(ltraces)>10000){
					ss <- seq(1,nrow(ltraces),by=nrow(ltraces)/10000)
					ltraces <-  ltraces[ss,]
				}
				colors <- rainbow(ncol(ltraces))
				ret <- list(retlist=retlist,ltraces=ltraces)
				## plot results
				pdf(pdf,onefile=TRUE)
				plot(as.numeric(rownames(ltraces)),ltraces[,1],type="l",xlab="iteration",ylab="Score",ylim=range(ltraces,na.rm=TRUE),col=colors[1],main="Score traces")
				if(ncol(ltraces)>1)
					sapply(2:ncol(ltraces), function(j,ltraces,colors) lines(as.numeric(rownames(ltraces)),ltraces[,j],col=colors[j]), ltraces=ltraces, colors=colors)
				## get the final network from all cores inferences
				#net <- retlist[[1]]$phi
				#net <- retlist[[1]]
				for(netnr in 2:cores) {
					ret2 <- retlist[[netnr]]
					plotdetailed(ret2$phi,stimuli=ret2$stimuli,weights=ret2$weights)
					#net2 <- retlist[[netnr]]
					#net$conf.act <- net$conf.act + net2$conf.act
					#net$conf.inh <- net$conf.inh + net2$conf.inh					
#					sel <- net$conf.act!=0 & net2$conf.act!=0
#					net$conf.act[sel] <- colMeans(rbind(net$conf.act[sel],net2$conf.act[sel]))
#					sel2 <- setdiff(1:length(net$phi),which(sel))
#					net$conf.act[sel] <- net$conf.act[sel2] + net2$conf.act[sel2]				
#					sel <- net$conf.inh!=0 & net2$conf.inh!=0
#					net$conf.inh[sel] <- colMeans(rbind(net$conf.inh[sel],net2$conf.inh[sel]))
#					sel2 <- setdiff(1:length(net$phi),which(sel))
#					net$conf.inh[sel] <- net$conf.inh[sel2] + net2$conf.inh[sel2]
				}
				#net$conf.act <- net$conf.act/cores
				#net$conf.inh <- net$conf.inh/cores
				#ret2 <- get.phi.final(net,th=th)
				#plotdetailed(ret2$phi,stimuli=ret2$stimuli,weights=ret2$weights)
				dev.off()
			} else {
				#phistart <- matrix(sample(c(0,1,2),n*n,replace=T), nrow=n, ncol=n, dimnames=list(phinames,phinames))
				if(is.null(phi))
					phistart <- matrix(0, nrow=n, ncol=n, dimnames=list(phinames,phinames))
				else
					phistart <- phi
				ret <- mcmc_ddepn(dat, phiorig=phiorig, phi=phistart, stimuli=stimuli,
						th=th, multicores=multicores, pdf=pdf, maxiterations=maxiterations,
						usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z,
						maxiter=maxiter,fanin=fanin, gam=gam, it=it, K=K)
			}
		}
	}
	ret
}

