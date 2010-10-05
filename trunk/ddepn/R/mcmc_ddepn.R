# TODO: Add comment
# 
# Author: benderc
###############################################################################
runmcmc <- function(x,dat,phiorig,phi,stimuli,th,multicores,outfile,maxiterations,
		usebics,cores,lambda,B,Z,samplelambda,hmmiterations,fanin,gam,it,K,burnin,
		priortype) {
	ret <- mcmc_ddepn(dat, phiorig=phiorig, phi=x$phi, stimuli=stimuli,
			th=th, multicores=multicores, outfile=x$outfile, maxiterations=maxiterations,
			usebics=usebics, cores=cores, lambda=lambda, B=B, Z=Z, samplelambda=samplelambda,
			hmmiterations=hmmiterations,fanin=fanin, gam=gam, it=it, K=K,
			burnin=burnin,priortype=priortype)
	ret
}

mcmc_ddepn <- function(dat, phiorig=NULL, phi=NULL, stimuli=NULL,
		th=0.8, multicores=FALSE, outfile=NULL, maxiterations=10000,
		usebics=FALSE, cores=2, lambda=NULL, B=NULL,Z=NULL,
		samplelambda=TRUE, hmmiterations=30, fanin=4,
		gam=NULL, it=NULL, K=NULL, burnin=1000,priortype="laplaceinhib") {
	if(!is.null(B))
		diag(B) <- 0
	if(!priortype %in% c("laplaceinhib","laplace","scalefree"))
		stop("Error, priortype must be one of 'laplaceinhib', 'laplace' or 'scalefree'.")
#	laplaceinhib <- laplace <- scalefree <- none <- FALSE
#	switch(priortype,
#			laplaceinhib = laplaceinhib<-TRUE,
#			laplace = laplace<-TRUE,
#			scalefree = scalefree<-TRUE,
#			none<-TRUE)
	#browser()
	#laplace <- !is.null(lambda) && !is.null(B) && !is.null(Z)
	#scalefree <- !is.null(gam) && !is.null(it) && !is.null(K)
	# initialise
	#dat[is.na(dat)] <- 0
	antibodies <- rownames(dat)
	tps <- unique(sapply(colnames(dat), function(x) strsplit(x,"_")[[1]][2]))
	reps <- table(sub("_[0-9].*$","",colnames(dat))) / length(tps)
	gammaposs <- propagate.effect.set(phi,stimuli)
#	longprop <- 1:max(length(tps),(nrow(phi)*100))
#	gammaposs <- propagate.effect.set(phi,longprop,stimuli,reps=reps)
#	gammaposs <- uniquegammaposs(gammaposs)
	# now get an initial gamma matrix
	gammax <- NULL
	for(sti in 1:length(stimuli)) {
		st <- stimuli[[sti]]
		indices <- grep(paste("^",paste(names(st),collapse="&"),"_",sep=""),colnames(gammaposs))
		gx <- replicatecolumns(gammaposs[,sort(sample(indices,length(tps),replace=TRUE))],reps[sti])
		gammax <- cbind(gammax, gx)
	}
	#gammax <- propagate.effect.set(phi,tps,stimuli,reps=reps)
	Ltmp <- likl(dat,gammax)
	Linit <- Ltmp$L
	thetax <- Ltmp$theta
	Linit[Linit==Inf] <- 0
	Linit[Linit==-Inf] <- 0
	Linit <- sum(Linit)
	aicinit <- get.aic(phi,Linit)
	bicinit <- get.bic(phi,Linit, length(dat))
	prinit <- prior(phi, lambda, B, Z, gam, it, K, priortype)
	#if(laplace || scalefree || laplaceinhib) {
	if(priortype %in% c("laplaceinhib","laplace","scalefree")) {
		posteriorinit <- Linit + prinit
		#posteriorinit <- posterior(phi, sum(Linit), lambda, B, Z, gam, it, K, priortype)
	} else {
		posteriorinit <- NULL
	}
	#}
	movetypes <- c("switchtype","delete","addactivation","addinhibition","revert") ## v1
	#movetypes <- c("add","delete","revert") ## v2
	
	bestmodel <- list(phi=phi,L=Linit,aic=aicinit,bic=bicinit,posterior=posteriorinit,dat=dat,
			theta=thetax, gamma=gammax, gammaposs=gammaposs, tps=tps, stimuli=stimuli, reps=reps,
			hmmiterations=hmmiterations, TSA=NULL, Tt=NULL, lastmove="addactivation", coords=c(1,1),
			lambda=lambda,B=B,Z=Z,pegm=1,pegmundo=1,nummoves=length(movetypes),fanin=fanin,
			gam=gam,it=it,K=K,phi.orig=phiorig, burnin=burnin,priortype=priortype,pr=prinit)#,
			#samplelambda=samplelambda)

	it <- 1
	stats <- matrix(0, nrow=maxiterations, ncol=15, dimnames=list(1:maxiterations, c("MAP", "tp","tn","fp","fn","sn","sp","lambda","acpt","lacpt","stmove","lratio","prratio","postratio","proposalratio")))
	freqa <- freqi <- eoccur <- bestmodel$phi
	freqa[freqa!=0] <- 0
	freqi[freqi!=0] <- 0
	eoccur[eoccur!=0] <- 0
	while(it < maxiterations) {
		cat("iteration ", it, " ")
		#if(laplace) {
		if(priortype=="laplaceinhib" || priortype=="laplace") {			
			if(samplelambda) {
				newlambda <- runif(1, bestmodel$lambda-1, bestmodel$lambda+1)
				newlambda <- min(max(0.01,newlambda),500)
			} else {
				newlambda <- bestmodel$lambda
			}
		#} else if(scalefree) {
		} else if(priortype=="scalefree") {
			#newgam <- runif(1, bestmodel$gam-1, bestmodel$gam+1)
			#newgam <- min(max(2,newgam),30) # gamma mustn't be smaller than 2
			newgam <- bestmodel$gam
		}
		#print(paste("+++++",newlambda, samplelambda))
		movetype <- sample(1:length(movetypes),1)
		if(all(bestmodel$phi==0))
			movetype <- sample(c(3,4),1) ## v1
			#movetype <- 1 # add, ## v2
		if(all(bestmodel$phi!=0))
			movetype <- sample(c(1,2,5),1)## v1
			#movetype <- sample(c(2,3),1) # delete, revert, ## v2
			#movetype <- sample(c(2,5),1)## v1
			
		st <- system.time(b1 <- mcmc_move(bestmodel, movetypes[movetype]))
		if(b1[[1]]$posterior==Inf || b1[[1]]$posterior==-Inf){
			print("Posterior of proposal is Inf. Please check.")
			browser()
		}
		if(priortype=="laplace" | priortype=="laplaceinhib") {
			ret <- mcmc_accept(bestmodel, b1, newlambda)
		} else if (priortype=="scalefree") {
			ret <- mcmc_accept(bestmodel, b1, newgam)
		}
		if(ret$bestproposal$posterior==Inf || ret$bestproposal$posterior==-Inf) {
			print("Posterior of accepted model is Inf. Please check.")
			browser()
		}
		liklihoodratio <- bestmodel$L / b1[[1]]$L
		priorratio <- bestmodel$pr / b1[[1]]$pr
		posteriorratio <- bestmodel$posterior / b1[[1]]$posterior
		proposalratio <- b1[[1]]$pegmundo / b1[[1]]$pegm
		bestmodel <- ret$bestproposal
	#	if(it>=burnin) {
			## count how often any edge occurred at a given position
			tmp <- bestmodel$phi
			tmp[bestmodel$phi==2] <- 0
			freqa <- freqa + tmp
			tmp <- bestmodel$phi
			tmp[bestmodel$phi==1] <- 0
			freqi <- freqi + (tmp/2)
			### get a 'final' network
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
		
			#lst <- get.phi.final(bestmodel,th=th)
			lst <- get.phi.final.mcmc(list(bestmodel), it, prob=.333, qu=.99999)[[1]]	
			if(!is.null(phiorig)) {
				comp <- compare.graphs.tc(phiorig=phiorig,phi=lst$phi)
			} else {
				comp <- rep(0,8)
				names(comp) <- c("tp","tn","fp","fn","sn","sp","prec","f1")
			}
#browser()
			if(priortype=="laplace" || priortype=="laplaceinhib") {
				replace <- as.matrix(unlist(c(bestmodel$posterior, comp[1:6], bestmodel$lambda, ret$acpt, ret$lacpt, st[3],liklihoodratio,priorratio,posteriorratio,proposalratio)))
			} else if(priortype=="scalefree") {
				replace <- as.matrix(unlist(c(bestmodel$posterior, comp[1:6], bestmodel$gam, ret$acpt, ret$lacpt, st[3],liklihoodratio,priorratio,posteriorratio,proposalratio)))
			}
			if(nrow(replace)!=1)
				replace <- t(replace)
			stats[it,] <- replace
	#	}		
		# some convergence statistic
		#SSW <- sd(stats[,"MAP"],na.rm=T)
		#SSB <- 
		#Rhat <- ((nrow(stats)-1)/nrow(stats)) * SSW
		#l <- maxiterations
		l <- it #nrow(stats)
		if(it%%250==1 && it >= burnin){
			if(!is.null(outfile))
				pdf(outfile)
			par(mfrow=c(3,3),mar=c(3,4,1,1),oma=c(1,1,1,1))
			plot(1:l, stats[1:l,"MAP"], type='l',ylab="",xlab="iteration",main="log posterior")
			perf <- mcmc_performance(bestmodel)
			hist(stats[1:l,"acpt"],breaks=100,main="acpt")
			hist(stats[1:l,"lacpt"],breaks=100,main="lacpt")
			hist(stats[1:l,"lambda"],breaks=100,main="lambda")
			if(is.null(phiorig)) {
				plot.new()
				text(0.5,0.5,labels="no origininal network given")
			} else {
				plotdetailed(phiorig,stimuli=bestmodel$stimuli,fontsize=25)
			}
			#if(it>burnin) {
				boxplot(as.data.frame(stats[(burnin:it),c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
						main=paste("avgSN: ", signif(median(stats[(burnin:it),"sn"]),digits=4), "avgSP: ", signif(median(stats[(burnin:it),"sp"]),digits=4)))
			#} else {
			#	boxplot(as.data.frame(stats[1:it,c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
			#			main=paste("BURNIN avgSN: ", signif(median(stats[,"sn"]),digits=4), "avgSP: ", signif(median(stats[,"sp"]),digits=4)))
			#}
			# partial autocorrelation function:
			R <- acf(stats[1:l,"MAP"])
			weights <- bestmodel$phi
			weights[bestmodel$phi==1] <- freqa[bestmodel$phi==1]
			weights[bestmodel$phi==2] <- freqi[bestmodel$phi==2]
			plotdetailed(bestmodel$phi,stimuli=bestmodel$stimuli,weights=bestmodel$weights,fontsize=15)
			if(!is.null(outfile))
				dev.off()
		}
		if(it%%1000==1 && !is.null(outfile)) {
			rdfile <- sub(".pdf$","_mcmcdata.RData",outfile)
			save(bestmodel,stats,freqa,freqi,it,file=rdfile)
		}
		it <- it + 1

	}
	bestmodel[["stats"]] <- stats
	bestmodel[["freqa"]] <- freqa
	bestmodel[["freqi"]] <- freqi
	gc(verbose=FALSE)
	bestmodel
}