# TODO: Add comment
# 
# Author: benderc
###############################################################################
runmcmc <- function(x,dat,phiorig,phi,stimuli,th,multicores,outfile,maxiterations,
		usebics,cores,lambda,B,Z,samplelambda,hmmiterations,fanin,gam,it,K,burnin,
		priortype,plotresults=TRUE,always_sample_sf=FALSE,scale_lik=FALSE, allow.stim.off=TRUE,
		debug=0, retobj=NULL) {
	ret <- mcmc_ddepn(dat, phiorig=phiorig, phi=x$phi, stimuli=stimuli,
			th=th, multicores=multicores, outfile=x$outfile, maxiterations=maxiterations,
			usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, samplelambda=samplelambda,
			hmmiterations=hmmiterations,fanin=fanin, gam=gam, it=it, K=K,
			burnin=burnin,priortype=priortype,plotresults=plotresults,always_sample_sf=always_sample_sf,scale_lik=scale_lik, allow.stim.off=allow.stim.off,
			debug=debug, retobj=x$retobj)
	ret
}

mcmc_ddepn <- function(dat, phiorig=NULL, phi=NULL, stimuli=NULL,
		th=0.8, multicores=FALSE, outfile=NULL, maxiterations=10000,
		usebics=FALSE, cores=2, lambda=NULL, B=NULL,Z=NULL,
		samplelambda=NULL, hmmiterations=30, fanin=4,
		gam=NULL, it=NULL, K=NULL, burnin=1000,priortype="laplaceinhib",plotresults=TRUE,
		always_sample_sf=FALSE,scale_lik=FALSE, allow.stim.off=TRUE,debug=0,retobj=NULL) {
	if(!is.null(outfile))
		outfile <- sub("\\.pdf","_stats.pdf", outfile)
	if(!is.null(B))
		diag(B) <- 0
	if(!priortype %in% c("laplaceinhib","laplace","scalefree","uniform"))
		stop("Error, for MCMC, usebics must be FALSE and priortype one out of 'laplaceinhib', 'laplace', 'scalefree' or 'uniform'.")
	antibodies <- rownames(dat)
	tps <- unique(sapply(colnames(dat), function(x) strsplit(x,"_")[[1]][2]))
	reps <- table(sub("_[0-9].*$","",colnames(dat))) / length(tps)
	gammaposs <- propagate.effect.set(phi,stimuli, allow.stim.off=allow.stim.off)
	# now get an initial gamma matrix
	if(!is.null(retobj)) {
		gammax <- retobj$gamma
	} else {
		gammax <- NULL
		for(sti in 1:length(stimuli)) {
			st <- stimuli[[sti]]
			indices <- grep(paste("^",paste(names(st),collapse="&"),"_",sep=""),colnames(gammaposs))
			gx <- replicatecolumns(gammaposs[,sort(sample(indices,length(tps),replace=TRUE))],reps[sti])
			gammax <- cbind(gammax, gx)
		}
	}
	Ltmp <- likl(dat,gammax,scale_lik)
	Linit <- Ltmp$L
	thetax <- Ltmp$theta
	Linit[Linit==Inf] <- 0
	Linit[Linit==-Inf] <- 0
	Linit <- sum(Linit,na.rm=TRUE)
	aicinit <- get.aic(phi,Linit)
	bicinit <- get.bic(phi,Linit, length(dat))
	prinit <- prior(phi, lambda, B, Z, gam, it, K, priortype)
	if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {
		posteriorinit <- Linit + prinit
	} else {
		posteriorinit <- NULL
	}
	## which types for the moves are needed?
	movetypes <- c("switchtype","delete","addactivation","addinhibition","revert","revswitch") ## v3
	if(!is.null(retobj)) {
		mu_run <- retobj$mu_run
		Qi <- retobj$Qi
		freqa <- retobj$freqa
		freqi <- retobj$freqi
		eoccur <- retobj$eoccur
		iter <- retobj$iter+1
		lastmove <- retobj$lastmove
	} else {
		mu_run <- 0
		Qi <- 0
		freqa <- freqi <- eoccur <- phi
		freqa[freqa!=0] <- 0
		freqi[freqi!=0] <- 0
		eoccur[eoccur!=0] <- 0
		iter <- 1
		lastmove <- "addactivation"
	}
	
	bestmodel <- list(phi=phi,phiorig=phiorig,L=Linit,aic=aicinit,bic=bicinit,posterior=posteriorinit,dat=dat,
			theta=thetax, gamma=gammax, gammaposs=gammaposs, tps=tps, stimuli=stimuli, reps=reps,
			hmmiterations=hmmiterations, lastmove=lastmove, coords=c(1,1),
			lambda=lambda,B=B,Z=Z,pegm=1,pegmundo=1,nummoves=length(movetypes),fanin=fanin,
			gam=gam,it=it,K=K,phi.orig=phiorig, burnin=burnin,priortype=priortype,pr=prinit
			,mu_run=mu_run,Qi=Qi,sd_run=NA,freqa=freqa,freqi=freqi,eoccur=eoccur,scalefac=0.005,scale_lik=scale_lik,
			allow.stim.off=allow.stim.off,iter=iter,samplelambda=samplelambda,always_sample_sf=always_sample_sf)
	## setup a matrix holding the statistics
	## TODO if thin==TRUE, this matrix
	## is of size maxiterations/x=10000, i.e. store every xth element
	if(!is.null(retobj)) {
		stats <- retobj$stats
		maxiterations <- iter + maxiterations -1 ## get new maximum iteration number
		stats <- rbind(stats, matrix(0, nrow=(maxiterations-iter+1), ncol=ncol(stats), dimnames=list(iter:maxiterations, colnames(stats))))
	} else {
		stats <- matrix(0, nrow=maxiterations, ncol=18,
			dimnames=list(1:maxiterations, c("MAP", "tp","tn","fp","fn","sn","sp",
							"lambda","acpt","lacpt","stmove",
							"lratio","prratio","postratio","proposalratio",
							"prior","liklihood","scalefac")))
	}
	##########################
	## MCMC chain
	##########################
	cls()
	if(debug==0) {
		txt <- paste("# Using inhibMCMC,",cores,"parallel chains #")
		pad <- paste(rep("#",nchar(txt)),collapse="")
		print(pad)
		print(txt)
		print(pad)
		pb <- txtProgressBar(min = 0, max = maxiterations, style = 3)
	}
	while(iter <= maxiterations) {
		if(debug==0) {
			setTxtProgressBar(pb, iter)
		}
		if(priortype %in% c("laplaceinhib", "laplace", "uniform")) {
			## if samplelambda tells to hold lambda "fixed" or to integrate, 
			## then don't change lambda
			#if(samplelambda=="fixed" | samplelambda=="integrate") {
			if(is.null(samplelambda)) {
				newlambda <- bestmodel$lambda
			} else { ## sample the lambda value, if a numeric value is given for samplelambda, use
				## 2*samplelambda as sampling interval around the current lambda
				## to sample a new lambda from, otherwise use 0.1 as interval size
				if(mode(samplelambda)=="numeric") {
					stepsize <- samplelambda
				} else {
					stepsize <- 0.1					
				}
				newlambda <- runif(1, bestmodel$lambda-stepsize, bestmodel$lambda+stepsize)
				newlambda <- min(max(1e-8,newlambda),100)
			}
		} else if(priortype=="scalefree") { ## keep gam fixed
			#newgam <- runif(1, bestmodel$gam-1, bestmodel$gam+1)
			#newgam <- min(max(2,newgam),30) # gamma mustn't be smaller than 2
			newgam <- bestmodel$gam
		}
		b1 <- movesproposed <- NULL
		while(is.null(b1)) {
			movetype <- sample(setdiff(movetypes,movesproposed),1)
			if(all(bestmodel$phi==0)) {
				movetype <- sample(setdiff(movetypes[c(3,4)],movesproposed),1) #v4
			}
			if(all(bestmodel$phi!=0)) {
				movetype <- sample(setdiff(movetypes[c(1,2,5,6)],movesproposed),1)#v4
			}
			movesproposed <- c(movesproposed, movetype)
			st <- system.time(b1 <- mcmc_move(bestmodel, movetype))##v4
		}
		if(b1[[1]]$posterior==Inf || b1[[1]]$posterior==-Inf){
			print("Posterior of proposal is Inf. Please check.")
			browser()
		}
		bestmodel[["iter"]] <- iter
		if(priortype %in% c("laplaceinhib","laplace","uniform")) {
			ret <- mcmc_accept(bestmodel, b1, newlambda, debug)
		} else if (priortype=="scalefree") {
			ret <- mcmc_accept(bestmodel, b1, newgam, debug)
		}
		if(ret$bestproposal$posterior==Inf || ret$bestproposal$posterior==-Inf) {
			print("Posterior of accepted model is Inf. Please check.")
			browser()
		}
		liklihoodratio <- bestmodel$L - b1[[1]]$L
		priorratio <- bestmodel$pr - b1[[1]]$pr
		posteriorratio <- bestmodel$posterior - b1[[1]]$posterior
		proposalratio <- b1[[1]]$pegmundo - b1[[1]]$pegm
		bestmodel <- ret$bestproposal
		## during burnin: find scalefac parameter to adjust the acceptance rate
		### -> perhaps the scalefactor should be calculated all the time,
		### in the profiles, there seems to be a kink in the likelihood profile when the 
		### sf is fixed after the burnin
		## to lie around 0.4. don't know if this is a reasonable level for acceptance rates,
		## suggested in Gelman 2003, chapter 11.10, recommended posterior simulation strategy
		## sample the scaling factor up to half of the burn-in, then hold fixed and let 
		## the sampler mix until the end of the burn-in
		if(iter<=(burnin/2) | always_sample_sf==TRUE) {
			## find scale factor that holds acpt around .4, see gelman 2003 for explanation
			if((posteriorratio+proposalratio)==0)
				if(iter==1 || all(stats[1:iter,"scalefac"]==Inf)) ## some fallback scalefactor
					scalefac <- 0.005
				else
					scalefac <- median(stats[1:iter,"scalefac"])
			else
				scalefac <- min(abs(log(0.4))/abs(posteriorratio+proposalratio),1) ## is kept smaller 1, log 0.4 yields sfs around .6
		} else {
			sf <- stats[1:(min(iter,burnin)),"scalefac"]
			scalefac <- max(1e-20,median(sf[sf!=Inf],na.rm=TRUE))
		}
		bestmodel$scalefac <- scalefac
		bestmodel[["iter"]] <- iter
		
		if(iter>burnin) {
			## count how often any edge occurred at a given position
			tmp <- bestmodel$phi
			tmp[bestmodel$phi==2] <- 0
			freqa <- freqa + tmp
			tmp <- bestmodel$phi
			tmp[bestmodel$phi==1] <- 0
			freqi <- freqi + (tmp/2)
			## get confidences for edges from the samplings
			eoccur <- freqi + freqa
			lst <- bestmodel
			conf.act <- freqa/eoccur
			conf.act[is.na(conf.act)] <- 0
			conf.inh <- freqi/eoccur
			conf.inh[is.na(conf.inh)] <- 0
			bestmodel[["freqa"]] <- freqa
			bestmodel[["freqi"]] <- freqi
			bestmodel[["conf.act"]] <- conf.act 
			bestmodel[["conf.inh"]] <- conf.inh
			bestmodel[["eoccur"]] <- eoccur
			bestmodel[["phi.orig"]] <- phiorig
			bestmodel[["burnin"]] <- burnin
			## update mean and standard deviations of theta parameters
			bth <- bestmodel$theta
			bth[is.na(bth)] <- 0
			mu_run_plus1 <- bestmodel$mu_run + 1/iter * (bth - bestmodel$mu_run)
			Qiplus1 <- bestmodel$Qi + (bth - bestmodel$mu_run) * (bth - mu_run_plus1)
			bestmodel[["mu_run"]] <- mu_run_plus1
			bestmodel[["Qi"]] <- Qiplus1
			bestmodel[["sd_run"]] <- sqrt((1/(iter-1)) * Qiplus1) 
		}
		## get an intermediate network from the samplings
		if(iter>burnin)
			lst <- get.phi.final(bestmodel,th=th) # set th around 0.8
		else
			lst <- bestmodel # if in burnin, just use whatever is there
		if(!is.null(phiorig) & iter > burnin) {
			comp <- compare.graphs.tc(phiorig=phiorig,phi=lst$phi)
		} else {
			comp <- rep(0,8)
			names(comp) <- c("tp","tn","fp","fn","sn","sp","prec","f1")
		}
		## save some statistics for this iteration
		if(priortype %in% c("laplaceinhib","laplace","uniform")) {
			replace <- as.matrix(unlist(c(lst$posterior, comp[1:6], lst$lambda, ret$acpt, ret$lacpt, st[3],liklihoodratio,priorratio,posteriorratio,proposalratio,lst$pr,lst$L,scalefac)))
		} else if(priortype=="scalefree") {
			replace <- as.matrix(unlist(c(lst$posterior, comp[1:6], lst$gam, ret$acpt, ret$lacpt, st[3],liklihoodratio,priorratio,posteriorratio,proposalratio,lst$pr,lst$L,scalefac)))
		}
		if(nrow(replace)!=1)
			replace <- t(replace)
		stats[iter,] <- replace
		# some convergence statistic -> should become gelmans Rhat at some time
		#SSW <- sd(stats[,"MAP"],na.rm=T)
		#SSB <- 
		#Rhat <- ((nrow(stats)-1)/nrow(stats)) * SSW
		if(iter%%250==1 && iter > burnin){
			if(plotresults) {
				if(!is.null(outfile))
					pdf(outfile,width=10,height=10)
				start <- burnin + 1
				layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = TRUE))
				## posterior
				plot(1:iter, stats[1:iter,"MAP"], type='l', ylab="", xlab="iteration", main="Posterior trace")
				abline(v=start,col="green")
				plot(1:iter, stats[1:iter,"postratio"], type='l', ylab="", xlab="iteration", main="Posterior ratios")
				abline(v=start,col="green")	
				## orig/inferred network
				if(is.null(phiorig)) {
					plot.new()
					text(0.5,0.5,labels="no origininal network given")
				} else {
					plotdetailed(phiorig,stimuli=lst$stimuli,fontsize=15,main="Original net")
				}
				## liklihood
				plot(1:iter, stats[1:iter,"liklihood"], type='l', ylab="", xlab="iteration", main="Liklihood trace")
				abline(v=start,col="green")
				plot(1:iter, stats[1:iter,"lratio"], type='l', ylab="", xlab="iteration", main="Liklihood ratios")
				abline(v=start,col="green")
				## inferred network
				plotdetailed(lst$phi,stimuli=lst$stimuli,weights=lst$weights,fontsize=15, main="Inferred net")	
				## prior
				plot(1:iter, stats[1:iter,"prior"], type='l', ylab="", xlab="iteration", main="Prior trace")
				abline(v=start,col="green")
				plot(1:iter, stats[1:iter,"prratio"], type='l', ylab="", xlab="iteration", main="Prior ratios")
				abline(v=start,col="green")
				## roc curve
				perf <- mcmc_performance(lst)
				# some more statistics that could be plotted
				### acceptance rate
				#stpl <- stats[1:iter,"acpt"]
				#hist(stpl[stpl!=1],breaks=100,main="Acceptance rates (only != 1)")
				### sn/sp plot
				#boxplot(as.data.frame(stats[((burnin+1):iter),c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
				#		main=paste("avgSN: ", signif(median(stats[(burnin:iter),"sn"]),digits=4), "avgSP: ", signif(median(stats[(burnin:iter),"sp"]),digits=4)))
				# partial autocorrelation function:
				#R <- acf(stats[1:iter,"MAP"])
				if(!is.null(outfile))
					dev.off()
			}
		}
		if(iter%%1000==1 && !is.null(outfile)) {
			rdfile <- sub(".pdf$","_mcmcdata.RData",outfile)
			save(bestmodel,stats,freqa,freqi,iter,file=rdfile)
		}
		iter <- iter + 1
	}
	print("done.")
	bestmodel[["stats"]] <- stats
	bestmodel[["freqa"]] <- freqa
	bestmodel[["freqi"]] <- freqi
	gc(verbose=FALSE)
	bestmodel
}