# Calling function to invoke the modelling framework
# Data matrix dat in the following format:
# 	rows: antibodies/nodes of network
#	cols: experiments, labeled by EXPLABEL_time
# lambda: laplace prior hyperparameter
# B     : prior weights matrix
# gam   : scalefree prior hyperparameter, also as exponent in
#         P(Phi|lambda,gam,B) = 1/2lambda * exp (-|Phi-B|^gam/lambda)
# it    : scalefree prior hyperparameter
# K     : scalefree prior hyperparameter
# fanin: maximum number of incoming edges, used for efficient computation of the 
#         normalization factor for the prior distribution
#
#
#  Priors: laplace: Froehlich 2007
#          scalefree: Kamimura and Shimodaira, A Scale-free Prior over Graph Structures for Bayesian Inference of Gene Networks
# Author: benderc
###############################################################################

ddepn <- function(dat, phiorig=NULL, phi=NULL, th=0.5, inference="netga", outfile=NULL,
                  multicores=FALSE, maxiterations=1000, p=500, q=0.3, m=0.8, P=NULL,
				  usebics=TRUE, cores=2, 
				  priortype="laplaceinhib",
				  lambda=NULL, B=NULL, samplelambda=TRUE,
				  hmmiterations=100, fanin=4,
				  gam=NULL,it=NULL,K=NULL,quantL=.5,quantBIC=.5,
				  debug=FALSE,burnin=1000, thin=FALSE) {
	# get the experiments, i.e. the stimuli/inhibitor combinations
	# works if format of dat is like:
	# colnames contain the experiments in form STIMULUS_time
	cols <- colnames(dat)
	tmp <- sapply(cols, function(x) strsplit(x,"_")[[1]])
	## check if number of time points is the same across all experiments
	if(length(unique(table(tmp[2,])))!=1) {
		stop("ERROR: Found differing number of time points across experiments.")
	}
	exps <- unique(tmp[1,])
	stims <- sapply(exps, function(x) strsplit(x,"&")[[1]])
	allstim <- unique(unlist(stims))
	stimuli <- list()
	for(i in 1:length(stims)) {
		el <- stims[[i]]
		# find the row in which the stimulus is in the data matrix
		x <- match(el,rownames(dat))
		# or define a number 
		if(any(is.na(x)))
			x <- match(el,allstim)
		names(x) <- el
		stimuli[[i]] <- x
	}
	# add the stimuli as dummy data rows, if they are missing
	stimm <- match(unique(names(unlist(stimuli))),rownames(dat))
	if(any(is.na(stimm))) {
		xx <- unlist(stimuli)
		xxmat <- unique(cbind(xx,names(xx)))
		toattach <- matrix(0.0,nrow=nrow(xxmat),ncol=ncol(dat),dimnames=list(xxmat[,2],colnames(dat)))
		## remove stimulus that is already there, to get the order right
		prune <- stimm[!is.na(stimm)]
		if(length(prune)>0)
			dat <- dat[-prune,]
		dat <- rbind(toattach,dat)
	}
	n <- nrow(dat)
	phinames <- rownames(dat)
	
	## check if seeds are provided properly
	if(!is.null(phi)) {
		# set up the seed population for netga
		if(inference=="netga") {
			if(class(phi)=="list" && length(phi)!=p) {
				stop(paste("Error: length of seed network list must be the same as p =",p))
			} else if(class(phi)=="matrix") {
				if(dim(phi)!=c(nrow(dat),nrow(dat))) {
					stop(paste("Error: dimension of seed network must be",nrow(dat),"x",nrow(dat),"."))
				}
				tmp <- vector("list",p)
				for(i in 1:p)
					tmp[[i]] <- phi
				phi <- tmp
				rm(tmp)
			} else if(is.null(phi)) {
				phi <- matrix(0, nrow=nrow(dat), ncol=ncol(dat), dimnames=list(rownames(dat),rownames(dat)))
				tmp <- vector("list",p)
				for(i in 1:p)
					tmp[[i]] <- phi
				phi <- tmp
				rm(tmp)
			} else {
				stop("Error: please provide either a list of or a single seed network in argument phi, or leave it as NULL.")
			}
			V <- rownames(dat)
			tps <- unique(sapply(colnames(dat), function(x) strsplit(x,"_")[[1]][2]))
			reps <- table(sub("_[0-9].*$","",colnames(dat))) / length(tps)
			X <- vector("list",p)
			for(i in 1:p) {
				X[[i]] <- phi[[i]]
			}
			if(multicores) {
				P <- mclapply(X, getfirstphi, dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K,priortype=priortype, mc.preschedule=FALSE,mc.cores=cores)		
			} else {
				P <- lapply(X, getfirstphi, dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K,priortype=priortype)
			}
		} else if(inference=="mcmc") {
			## if a list of networks is given, create a named list that is used for mcmc_ddepn
			if(class(phi)=="list"){
				if(length(phi)!=cores)
					stop(paste("Error: length of seed network list must be the same as cores =",cores))
				P <- vector("list",cores)
				for(i in 1:cores) {
					P[[i]] <- list(phi=phi[[i]])
				}
			} else if(class(phi)=="matrix") {
				## if a matrix is given, create a list of cores copies
				## for cores independent mcmc runs.
				if(dim(phi)!=c(nrow(dat),nrow(dat))) {
					stop(paste("Error: dimension of seed network must be",nrow(dat),"x",nrow(dat),"."))
				}
				if(multicores==TRUE) {
					## create list of start nets, copy the 
					## one start network cores times
					P <- vector("list",cores)
					for(i in 1:cores)
						P[[i]] <- list(phi=phi)
					rm(tmp)
				}
			}
		}
	}
	## if BIC score should be used, don't use a prior
	if(usebics)
		priortype <- "none"

	## preparing the priors
	Z <- NULL
	if(priortype %in% c("laplace", "laplaceinhib") && !usebics) {
		if(is.null(lambda) | is.null(B))
			stop("Please specify arguments lambda and B for use of laplaceinhib or laplace prior.")
		## get normalisation factor for the networks
		##print("Computing prior normalisation factor...")
		##Z <- zlambda(B, lambda)
		##print("done.")
	} else if(priortype %in% c("scalefree") && !usebics) {
		if(is.null(gam) | is.null(it) | is.null(K))
			stop("Please specify arguments gam, it and K for use of scalefree prior.")
	}
	#else if(priortype %in% c("none")) {
	#	stop("Error in function arguments. Please specifiy either lambda/gamma for laplace prior, gam/it/K for scalefree prior or none if no prior distribution should be used.")
	#}
	## if GA should be used
	if(inference=="netga") {
		if(!is.null(outfile)) {
			bname <- basename(outfile)
			dname <- dirname(outfile)
			scorefile <- paste(dname,paste("score",sub("\\.pdf","",bname),".pdf",sep=""),sep="/")
		} else {
			scorefile <- NULL
		}
        stime <- system.time(retnetga <- netga(dat,stimuli,P=P,maxiterations=maxiterations,
						p=p,q=q,m=m,multicores=multicores,usebics=usebics,
						cores=cores,lambda=lambda,B=B,Z=Z,hmmiterations=hmmiterations,
						scorefile=scorefile,fanin=fanin,
						gam=gam,it=it,K=K,quantL=quantL,quantBIC=quantBIC,priortype=priortype))
		P <- retnetga$P
		scorestats <- retnetga$scorestats
		rm(retnetga)
		phi.activation.count <- phi.inhibition.count <- weights.tc <- matrix(0,nrow=nrow(dat),ncol=nrow(dat),dimnames=list(rownames(dat),rownames(dat))) 
		# now compare all graphs to the original
        result <- NULL
		resultep <- NULL
        for(i in 1:length(P)) {
			if(!is.null(phiorig))
            	result <- rbind(result, cbind(compare.graphs.tc(phiorig, P[[i]]$phi,ignore.type=FALSE),t(as.matrix(stime))))
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
					weights.tc=weights.tc, stats=result, conf.act=conf.act,conf.inh=conf.inh,
					stimuli=stimuli) #, resultep=resultep)
        ret <- get.phi.final(ret,th)
        plotrepresult(ret,outfile)
        ret[["P"]] <- P
		ret[["scorestats"]] <- scorestats
	} else {
		if(inference=="mcmc") {
			if(multicores) {
				## start an mcmc run for each core
				## list of start networks: if not given, create cores empty
				## networks
				if(is.null(P)) {
					P <- list()
					phistart <- matrix(0, nrow=n, ncol=n, dimnames=list(phinames,phinames))
					for(cr in 1:cores) {
						P[[cr]] <- list(phi=phistart)
					}
				}
				## prepare a temp object for the runmcmc call
				X <- list()
				for(cr in 1:cores) {
					if(!is.null(outfile))
						filename <- paste(sub("\\.pdf","",outfile),"_",cr,".pdf",sep="")
					else
						filename <- NULL
					X[[cr]] <- list(phi=P[[cr]]$phi, outfile=filename)
				}
				retlist <- mclapply(X, runmcmc, dat=dat, phiorig=phiorig, phi=NULL, stimuli=stimuli,
						th=th, multicores=multicores, outfile=outfile, maxiterations=maxiterations,
						usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, samplelambda=samplelambda,
						hmmiterations=hmmiterations,fanin=fanin, gam=gam, it=it, K=K, burnin=burnin,
						priortype=priortype,
						mc.preschedule=FALSE,mc.cores=cores)
				### experimental
				if(debug)
					browser()
				#retlist <- get.phi.final.mcmc(retlist, maxiterations, prob=.3333, qu=.99999)
				##### end experimental
			} else {
				#phistart <- matrix(sample(c(0,1,2),n*n,replace=T), nrow=n, ncol=n, dimnames=list(phinames,phinames))
				if(is.null(phi)) {
					## create empty start net, if phi==NULL
					phistart <- matrix(0, nrow=n, ncol=n, dimnames=list(phinames,phinames))
				} else {
					## otherwise, take the first net from phi
					if(class(phi)=="list")
						phistart <- phi[[1]]
					else
						phistart <- phi
				}
				ret <- mcmc_ddepn(dat, phiorig=phiorig, phi=phistart, stimuli=stimuli,
						th=th, multicores=multicores, outfile=outfile, maxiterations=maxiterations,
						usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, samplelambda=samplelambda,
						hmmiterations=hmmiterations,fanin=fanin, gam=gam, it=it, K=K, burnin=burnin,
						priortype=priortype)
				## convert the return value to a list, to keep it in the same format
				## as for the multicore==TRUE case
				retlist <- list()
				retlist[[1]] <- ret
			}
			## get the likelihood traces
			ltraces <- as.matrix(sapply(retlist,function(x) x$stats[,"MAP"]))
			## thinning: make sure ltraces has less than 10000 rows
			## TODO: change such that only 10000 rows are recorded (i.e. in mcmc_ddepn),
			## otherwise it doesn't really make sense
			if(nrow(ltraces)>10000 && thin == TRUE){
				ss <- seq(1,nrow(ltraces),by=nrow(ltraces)/10000)
				ltraces <-  as.matrix(ltraces[ss,])
			}
			colors <- rainbow(ncol(ltraces))
			ret <- list(samplings=retlist,ltraces=ltraces)
			## plot results
			if(!is.null(outfile)) {
				pdf(outfile,onefile=TRUE)
			} else {
				x11()
				par(mfrow=c(1,2))
			}
			plot(as.numeric(rownames(ltraces)),ltraces[,1],type="l",xlab="iteration",ylab="Score",ylim=range(ltraces,na.rm=TRUE),col=colors[1],main="Score traces")
			if(ncol(ltraces)>1)
				sapply(2:ncol(ltraces), function(j,ltraces,colors) lines(as.numeric(rownames(ltraces)),ltraces[,j],col=colors[j]), ltraces=ltraces, colors=colors)
			## get the final network from all cores inferences and plot
			for(netnr in 1:length(retlist)) {
				ret2 <- retlist[[netnr]]
				plotdetailed(ret2$phi,stimuli=ret2$stimuli,weights=ret2$weights)
			}
			if(!is.null(outfile))
				dev.off()
		}
	}
	ret
}

