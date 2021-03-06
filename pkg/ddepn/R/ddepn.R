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

ddepn <- function(dat, phiorig=NULL, phi=NULL, th=0.8, inference="netga", outfile=NULL,
                  multicores=FALSE, maxiterations=1000, p=500, q=0.3, m=0.8, P=NULL,
				  usebics=TRUE, cores=1, 
				  priortype="laplaceinhib",
				  lambda=NULL, B=NULL, samplelambda=NULL,
				  hmmiterations=100, fanin=4,
				  gam=NULL,it=NULL,K=NULL,quantL=.5,quantBIC=.5,
				  debug=0,burnin=500, thin=FALSE, plotresults=TRUE,
				  always_sample_sf=FALSE,scale_lik=FALSE,allow.stim.off=FALSE,
				  implementation="C") {
	## deal with missing replicates
	## make pad of NA columns for missing replicates
	dat <- pad_data(dat)
	## order the new matrix
	dat <- order_experiments(dat)		  
	# get the experiments, i.e. the stimuli/inhibitor combinations
	# works if format of dat is like:
	# colnames contain the experiments in form STIMULUS_time
	cols <- colnames(dat)
	tmp <- sapply(cols, function(x) strsplit(x,"_")[[1]])
	## check if number of time points is the same across all experiments
	## should be possible now
#	if(length(unique(table(tmp[2,])))!=1) {
#		stop("ERROR: Found differing number of time points across experiments.")
#	}
	## check if samplelambda is either numeric or NULL
	if(!(mode(samplelambda)=="numeric" | is.null(samplelambda))) {
		stop("ERROR: Please specify argument samplelambda either as numeric value or NULL.")	
	}
	## check if samplelambda is either numeric or NULL
	if(!is.null(lambda)) {
		if(!(mode(lambda)=="numeric" | is.na(lambda))) {
			stop("ERROR: Please specify argument lambda either as numeric value, NULL or NA.")	
		}
	}
	tmp2 <- addstimuli(dat)
	dat <- tmp2$dat
	stimuli <- tmp2$stimuli
	rm(tmp2)
	N <- nrow(dat)
	NC <- ncol(dat)
	V <- rownames(dat)
	
	## check if seeds are provided properly
	if(!is.null(phi)) {
		# set up the seed population for netga
		if(inference=="netga") {
			if(class(phi)=="list" && length(phi)!=p) {
				stop(paste("Error: length of seed network list must be the same as p =",p))
			} else if(class(phi)=="matrix") {
				if(dim(phi)!=c(N,N)) {
					stop(paste("Error: dimension of seed network must be",N,"x",N,"."))
				}
				tmp <- vector("list",p)
				for(i in 1:p)
					tmp[[i]] <- phi
				phi <- tmp
				rm(tmp)
			} else if(is.null(phi)) {
				phi <- matrix(0, nrow=N, ncol=NC, dimnames=list(rownames(dat),rownames(dat)))
				tmp <- vector("list",p)
				for(i in 1:p)
					tmp[[i]] <- phi
				phi <- tmp
				rm(tmp)
			} else {
				stop("Error: please provide either a list of or a single seed network in argument phi, or leave it as NULL.")
			}
			#V <- rownames(dat)
			ordstim <- sapply(stimuli, function(x) paste(names(x),collapse="&"))
			tmp <- tapply(colnames(dat), gsub("_.*$","",colnames(dat)), get_reps_tps)
			tmp <- tmp[ordstim]
			tps <- lapply(tmp, function(x) x$tps)
			reps <- unique(sapply(tmp, function(x) as.numeric(x$reps)))
			if(length(reps)>1) {
				stop("Number of replicates differ in the experiments. Please add columns containing NAs to obtain equal number of replicates.")
			}
			X <- vector("list",p)
			for(i in 1:p) {
				X[[i]] <- phi[[i]]
			}
			if(multicores) {
				P <- mclapply(X, getfirstphi, dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K,priortype=priortype,scale_lik=scale_lik, allow.stim.off=allow.stim.off, implementation=implementation, mc.preschedule=FALSE,mc.cores=cores)		
			} else {
				P <- lapply(X, getfirstphi, dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K,priortype=priortype,scale_lik=scale_lik, allow.stim.off=allow.stim.off, implementation=implementation)
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
				if(dim(phi)!=c(N,N)) {
					stop(paste("Error: dimension of seed network must be",N,"x",N,"."))
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
		## make sure the prior matrix only contains 0 and 1 entries when prior is 'laplace'
		if(priortype=="laplace")
			B <- detailed.to.simple.regulations(B)
		## get normalisation factor for the networks
		##print("Computing prior normalisation factor...")
		##Z <- zlambda(B, lambda)
		##print("done.")
	} else if(priortype %in% c("scalefree") && !usebics) {
		if(is.null(gam) | is.null(it) | is.null(K))
			stop("Please specify arguments gam, it and K for use of scalefree prior.")
	} else if(priortype %in% c("none","uniform")) {
		gam <- NULL
		it <- NULL
		K <- NULL
		B <- NULL
		lambda <- 0
	}
	## if integration should be performed over lambda, set
	## to NA, which is recognized in the prior function
	## and integration is performed there.
	#if(samplelambda=="integrate") {
	#	lambda <- NA
	#}
	## assign global variable IMPLEMENTATION
    ##envDDEPN <- new.env()
    ##assign("IMPLEMENTATION", implementation, envDDEPN)
	#assign("IMPLEMENTATION", implementation, .GlobalEnv)
	
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
						gam=gam,it=it,K=K,quantL=quantL,quantBIC=quantBIC,priortype=priortype,
						plotresults=plotresults,scale_lik=scale_lik, allow.stim.off=allow.stim.off,
						debug=debug, implementation=implementation))
		P <- retnetga$P
		scorestats <- retnetga$scorestats
		rm(retnetga)
		phi.activation.count <- phi.inhibition.count <- weights.tc <- matrix(0,nrow=N,ncol=N,dimnames=list(rownames(dat),rownames(dat))) 
		# now compare all graphs to the original
        result <- NULL
		resultep <- NULL
        for(i in 1:length(P)) {
			if(!is.null(phiorig))
            	result <- rbind(result, cbind(compareGraphs(phiorig, P[[i]]$phi,ignore.type=FALSE),t(as.matrix(stime))))
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
					stimuli=stimuli,
					quantBIC=quantBIC,quantL=quantL,q=q,m=m,usebics=usebics, implementation=implementation)
					#p=p,q=q,m=m,multicores=multicores,usebics=usebics,
					#cores=cores,lambda=lambda,B=B,Z=Z,hmmiterations=hmmiterations,
					#fanin=fanin, gam=gam,it=it,K=K,quantL=quantL,quantBIC=quantBIC,priortype=priortype,
					#scale_lik=scale_lik, allow.stim.off=allow.stim.off,debug=debug)
        ret <- get.phi.final(ret,th)
		if(plotresults)
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
					phistart <- matrix(0, nrow=N, ncol=N, dimnames=list(V,V))
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
						priortype=priortype, plotresults=plotresults,always_sample_sf=always_sample_sf,scale_lik=scale_lik,
						allow.stim.off=allow.stim.off,debug=debug, implementation=implementation,
						mc.preschedule=FALSE,mc.cores=cores)
				#retlist <- get.phi.final.mcmc(retlist, maxiterations, prob=.3333, qu=.99999)
				##### end experimental
			} else {
				#phistart <- matrix(sample(c(0,1,2),n*n,replace=T), nrow=n, ncol=n, dimnames=list(phinames,phinames))
				if(is.null(phi)) {
					## create empty start net, if phi==NULL
					phistart <- matrix(0, nrow=N, ncol=N, dimnames=list(V,V))
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
						priortype=priortype, plotresults=plotresults,always_sample_sf=always_sample_sf,scale_lik=scale_lik,
						allow.stim.off=allow.stim.off,debug=debug, implementation=implementation)
				## convert the return value to a list, to keep it in the same format
				## as for the multicore==TRUE case
				retlist <- list()
				retlist[[1]] <- ret
			}
			## get the likelihood traces
			ltraces <- as.matrix(sapply(retlist,function(x) x$stats[,"MAP"]))
			## thinning: make sure ltraces has less than 10000 rows
			## TODO: change such that only 10000 rows are recorded (i.e. in mcmc_ddepn),
			## otherwise it doesn't help to save memory for the trace storage
			if(nrow(ltraces)>10000 && thin == TRUE){
				ss <- seq(1,nrow(ltraces),by=nrow(ltraces)/10000)
				ltraces <-  as.matrix(ltraces[ss,])
			}
			colors <- rainbow(ncol(ltraces))
			ret <- list(samplings=retlist,ltraces=ltraces)
			## plot results
			if(!is.null(outfile)) {
				pdf(outfile,onefile=TRUE)
			}
			if(plotresults) {
				#if(is.null(outfile))
				#	x11()
				## find a suitable plot layout with little empty segments
				nplots <- length(retlist)+1
				mfdim1 <- mfdim2 <- ceiling(sqrt(nplots))
				if((mfdim1*mfdim2 - nplots) > mfdim1) {
					mfdim2 <- mfdim2 - 1
				} 
				par(mfrow=c(mfdim1,mfdim2))
				plot(as.numeric(rownames(ltraces)),ltraces[,1],type="l",xlab="iteration",ylab="Score",ylim=range(ltraces,na.rm=TRUE),col=colors[1],main="Score traces")
				if(ncol(ltraces)>1)
					sapply(2:ncol(ltraces), function(j,ltraces,colors) lines(as.numeric(rownames(ltraces)),ltraces[,j],col=colors[j]), ltraces=ltraces, colors=colors)
				## get the final network from all cores inferences and plot
				for(netnr in 1:length(retlist)) {
					ret2 <- retlist[[netnr]]
					plotdetailed(ret2$phi,stimuli=ret2$stimuli,weights=ret2$weights, main=paste("MCMC run", netnr))
				}
			}
			if(!is.null(outfile))
				dev.off()
		}
	}
	ret
}

## uses a returned object from netga or inhibMCMC and resumes the sampling/optimisation
resume_ddepn <- function(ret,maxiterations=10000,outfile=NULL,th=0.8,plotresults=TRUE,debug=0,cores=NULL, implementation="C", thin=FALSE) {
	#envDDEPN <- new.env()
    #assign("IMPLEMENTATION", implementation, envDDEPN)
	## close all x11 connections
	graphics.off()
	if(is.null(ret$samplings)) {
		if(is.null(cores))
			cores <- 1
		mc <- cores>1
		iter <- nrow(ret$scorestats)
		maxiterations <- iter + maxiterations
		## restore inference state
		p <- length(ret$P)
		q <- ret$q
		m <- ret$m
		quantBIC <- ret$quantBIC
		quantL <- ret$quantL
		usebics <- ret$usebics
		pp <- ret$P[[1]]
		lambda <- pp$lambda
		B <- pp$B
		Z <- pp$Z
		hmmiterations <- pp$hmmiterations
		fanin <- pp$fanin
		gam <- pp$gam
		it <- pp$it
		K <- pp$K
		priortype <- pp$priortype
		scale_lik <- pp$scale_lik
		allow.stim.off <- pp$allow.stim.off
		N <- nrow(pp$phi)
		dat <- ret$dat
		stimuli <- ret$stimuli
		phiorig <- ret$phi.orig
		P <- ret$P
		if(!is.null(outfile)) {
			bname <- basename(outfile)
			dname <- dirname(outfile)
			scorefile <- paste(dname,paste("score",sub("\\.pdf","",bname),".pdf",sep=""),sep="/")
		} else {
			scorefile <- NULL
		}
		stime <- system.time(retnetga <- netga(dat,stimuli,P=P,maxiterations=maxiterations,
						p=length(ret$P),q=q,m=m,multicores=mc,usebics=usebics,
						cores=cores,lambda=lambda,B=B,Z=Z,hmmiterations=hmmiterations,
						scorefile=scorefile,fanin=fanin,
						gam=gam,it=it,K=K,quantL=quantL,quantBIC=quantBIC,priortype=priortype,
						plotresults=plotresults,scale_lik=scale_lik, allow.stim.off=allow.stim.off,
						debug=debug,retobj=ret, implementation=implementation))
		P <- retnetga$P
		scorestats <- retnetga$scorestats
		rm(retnetga)
		phi.activation.count <- phi.inhibition.count <- weights.tc <- matrix(0,nrow=N,ncol=N,dimnames=list(rownames(dat),rownames(dat))) 
		# now compare all graphs to the original
		result <- NULL
		resultep <- NULL
		for(i in 1:length(P)) {
			if(!is.null(phiorig))
				result <- rbind(result, cbind(compareGraphs(phiorig, P[[i]]$phi,ignore.type=FALSE),t(as.matrix(stime))))
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
				stimuli=stimuli,
				quantBIC=quantBIC,quantL=quantL,q=q,m=m,usebics=usebics, implementation=implementation)
		ret <- get.phi.final(ret,th)
		if(plotresults)
			plotrepresult(ret,outfile)
		ret[["P"]] <- P
		ret[["scorestats"]] <- scorestats
	} else {
		if(!is.null(cores)) {
			print("NOTE: Argument cores is derived from the previous run. Omitting the given function argument cores.")
			Sys.sleep(3)
		}
		## find out if multicores have to be used
		cores <- ncol(ret$ltraces)
		mc <- cores>1
		## get the actual networks and the output file list
		X <- list()
		for(cr in 1:cores) {
			if(!is.null(outfile))
				filename <- paste(sub("\\.pdf","",outfile),"_",cr,".pdf",sep="")
			else
				filename <- NULL
			retobj <- ret$samplings[[cr]]
			X[[cr]] <- list(phi=retobj$phi, outfile=filename, retobj=retobj)
		}
		rs <- ret$samplings[[1]]
		usebics <- rs$posterior>0
		## phi is not used in runmcmc
		## outfile is not used in runmcmc
		## samplelambda in bestmodel! DONE
		## it in bestmodel ist für scalefree prior, iter als neues element für die iteration! DONE
		## always_sample_sf in bestmodel! DONE
		if(mc) {
			retlist <- mclapply(X, runmcmc, dat=rs$dat, phiorig=rs$phi.orig, phi=NULL, stimuli=rs$stimuli,
						th=th, multicores=mc, outfile=NULL, maxiterations=maxiterations,
						usebics=usebics, cores=cores, lambda=rs$lambda, B=rs$B, Z=rs$Z, samplelambda=rs$samplelambda,
						hmmiterations=rs$hmmiterations,fanin=rs$fanin, gam=rs$gam, it=rs$it, K=rs$K, burnin=rs$burnin,
						priortype=rs$priortype, plotresults=plotresults,always_sample_sf=rs$always_sample_sf,scale_lik=rs$scale_lik,
						allow.stim.off=rs$allow.stim.off,debug=debug, implementation=implementation,
						mc.preschedule=FALSE,mc.cores=cores)
		} else {
			ret <- runmcmc(X[[1]],rs$dat,rs$phi.orig,phi=NULL,rs$stimuli,th,mc,outfile=NULL,maxiterations,
						usebics,cores,rs$lambda,rs$B,Z=NULL,samplelambda=rs$samplelambda,rs$hmmiterations,rs$fanin,rs$gam,rs$it,rs$K,burnin=rs$burnin,
						rs$priortype,plotresults=plotresults,always_sample_sf=rs$always_sample_sf,rs$scale_lik,rs$allow.stim.off,debug=debug, 
						implementation=implementation)
			retlist <- list()
			retlist[[1]] <- ret	
		}
		## get the likelihood traces
		ltraces <- as.matrix(sapply(retlist,function(x) x$stats[,"MAP"]))
		## thinning: make sure ltraces has less than 10000 rows
		## TODO: change such that only 10000 rows are recorded (i.e. in mcmc_ddepn),
		## otherwise it doesn't help to save memory for the trace storage
		if(nrow(ltraces)>10000 && thin == TRUE){
			ss <- seq(1,nrow(ltraces),by=nrow(ltraces)/10000)
			ltraces <-  as.matrix(ltraces[ss,])
		}
		colors <- rainbow(ncol(ltraces))
		ret <- list(samplings=retlist,ltraces=ltraces)
		## plot results
		if(!is.null(outfile)) {
			pdf(outfile,onefile=TRUE)
		}
		if(plotresults) {
			plot(as.numeric(rownames(ltraces)),ltraces[,1],type="l",xlab="iteration",ylab="Score",ylim=range(ltraces,na.rm=TRUE),col=colors[1],main="Score traces")
			if(ncol(ltraces)>1)
				sapply(2:ncol(ltraces), function(j,ltraces,colors) lines(as.numeric(rownames(ltraces)),ltraces[,j],col=colors[j]), ltraces=ltraces, colors=colors)
			## get the final network from all cores inferences and plot
			for(netnr in 1:length(retlist)) {
				ret2 <- retlist[[netnr]]
				plotdetailed(ret2$phi,stimuli=ret2$stimuli,weights=ret2$weights)
			}
		}
		if(!is.null(outfile))
			dev.off()
	}
	return(ret)
}

