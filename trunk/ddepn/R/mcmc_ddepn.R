# TODO: Add comment
# 
# Author: benderc
###############################################################################


mcmc_ddepn <- function(dat, phiorig=NULL, phi=NULL, stimuli=NULL,
		th=0.8, multicores=FALSE, pdf=NULL, maxiterations=10000,
		usebics=FALSE, cores=2, lambda=NULL, B=NULL,Z=NULL,maxiter=30,fanin=4) {
	diag(B) <- 0
	# initialise
	#dat[is.na(dat)] <- 0
	antibodies <- rownames(dat)
	tps <- unique(sapply(colnames(dat), function(x) strsplit(x,"_")[[1]][2]))
	#reps <- ((ncol(dat)/length(tps))/length(stimuli))
	reps <- table(sub("_[0-9].*$","",colnames(dat))) / length(tps)
	longprop <- 1:max(length(tps),(nrow(phi)*100))
	gammaposs <- propagate.effect.set(phi,longprop,stimuli,reps=reps)
	gammaposs <- uniquegammaposs(gammaposs)
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
	if(is.null(lambda)) {
		posteriorinit <- NULL
	} else {
		posteriorinit <- posterior(phi, sum(Linit), lambda, B, Z)
	}
	#}
	#movetypes <- c("switchtype","delete","addactivation","addinhibition","revert") ## v1
	movetypes <- c("add","delete","revert") ## v2
	
	bestmodel <- list(phi=phi,L=Linit,aic=aicinit,bic=bicinit,posterior=posteriorinit,dat=dat,
			theta=thetax, gamma=gammax, gammaposs=gammaposs, tps=tps, stimuli=stimuli, reps=reps,
			maxiter=maxiter, TSA=NULL, Tt=NULL, lastmove="addactivation", coords=c(1,1), lambda=lambda,
			B=B,Z=Z,pegm=1,pegmundo=1,nummoves=length(movetypes),fanin=fanin)

	it <- 1
	stats <- matrix(0, nrow=maxiterations, ncol=11, dimnames=list(1:maxiterations, c("MAP", "tp","tn","fp","fn","sn","sp","lambda","acpt","lacpt","stmove")))
	freqa <- freqi <- eoccur <- bestmodel$phi
	freqa[freqa!=0] <- 0
	freqi[freqi!=0] <- 0
	eoccur[eoccur!=0] <- 0
	while(it < maxiterations) {
		cat("++ ", it, " ")
		newlambda <- runif(1, bestmodel$lambda-1, bestmodel$lambda+1)
		newlambda <- min(max(0.01,newlambda),30)
		movetype <- sample(1:length(movetypes),1)
		if(all(bestmodel$phi==0))
			movetype <- 1 # add, ## v2
			#movetype <- sample(c(3,4),1) ## v1
		if(all(bestmodel$phi!=0))
			movetype <- sample(c(2,3),1) # delete, revert, ## v2
			#movetype <- sample(c(2,5),1)## v1
			#movetype <- sample(c(1,2,5),1)## v1
		
		st <- system.time(b1 <- mcmc_move(bestmodel, movetypes[movetype]))
		ret <- mcmc_accept(bestmodel, b1, newlambda)
		bestmodel <- ret$bestproposal
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
		bestmodel[["conf.act"]] <- conf.act 
		bestmodel[["conf.inh"]] <- conf.inh
		bestmodel[["eoccur"]] <- eoccur
		bestmodel[["phiorig"]] <- phiorig
		lst <- get.phi.final(bestmodel,th=0.5)
		
		if(!is.null(phiorig)) {
			comp <- compare.graphs.tc(O=phiorig,M=lst$phi)
			stats[it,] <- as.matrix(cbind(bestmodel$posterior, comp[1:6], bestmodel$lambda, ret$acpt, ret$lacpt, st[3]))
		}
		# some convergence statistic
		#SSW <- sd(stats[,"MAP"],na.rm=T)
		#SSB <- 
		#Rhat <- ((nrow(stats)-1)/nrow(stats)) * SSW
		#l <- maxiterations
		l <- it #nrow(stats)
		if(it%%100==1 && it>1){
			if(!is.null(pdf))
				pdf(pdf)
			par(mfrow=c(3,3),mar=c(3,4,1,1),oma=c(1,1,1,1))
			plot(1:l, stats[1:l,"MAP"], type='l',ylab="",xlab="iteration",main="log posterior")
			perf <- mcmc_performance(bestmodel)
			#plot(1:l, stats[1:l,"sn"], type='l',ylim=c(0,1),ylab="",xlab="iteration",col="blue",main="SN/SP plot")
			#lines(1:l, stats[1:l,"sp"], lty=2,col="violet")
			#legend("topleft", c("sn","sp"), lty=c(1,2),col=c("blue","violet")) #,"red","orange"))
			hist(stats[1:l,"acpt"],breaks=100,main="acpt")
			hist(stats[1:l,"lacpt"],breaks=100,main="lacpt")
			hist(stats[1:l,"lambda"],breaks=100,main="lambda")
			#plot(1:l, stats[1:l,"acpt"], lty=1, col="#222222",main="acpt",pch=".")
			#plot(1:l, stats[1:l,"lacpt"], lty=1, col="#222222",main="lacpt",pch=".")
			#plot(1:l, stats[1:l,"lambda"], type='l',ylim=c(0,30),ylab="lambda",xlab="iteration")			
			#boxplot(as.data.frame(stats[,"stmove"]),ylab="time")
			plotdetailed(phiorig,stimuli=bestmodel$stimuli,fontsize=15)
			if(it>1000) {
				boxplot(as.data.frame(stats[(1000:it),c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
						main=paste("avgSN: ", signif(median(stats[(1000:it),"sn"]),digits=4), "avgSP: ", signif(median(stats[(1000:it),"sp"]),digits=4)))
			} else {
				boxplot(as.data.frame(stats[1:it,c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
						main=paste("BURNIN avgSN: ", signif(median(stats[,"sn"]),digits=4), "avgSP: ", signif(median(stats[,"sp"]),digits=4)))
			}
			# partial autocorrelation function:
			R <- acf(stats[1:l,"MAP"])
			weights <- bestmodel$phi
			weights[bestmodel$phi==1] <- freqa[bestmodel$phi==1]
			weights[bestmodel$phi==2] <- freqi[bestmodel$phi==2]
			plotdetailed(bestmodel$phi,stimuli=bestmodel$stimuli,weights=bestmodel$weights,fontsize=15)
			if(!is.null(pdf))
				dev.off()
		}
		it <- it + 1
		if(it%%500==0) {
			if(!is.null(pdf)) {
				save.image(file="MCMCimage.RData")	
			}
		}
	}
	bestmodel[["stats"]] <- stats
	bestmodel[["freqa"]] <- freqa
	bestmodel[["freqi"]] <- freqi
	bestmodel
}