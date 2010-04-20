# TODO: Add comment
# 
# Author: benderc
###############################################################################


mcmc_ddepn <- function(dat, phiorig=NULL, phi=NULL, stimuli=NULL,
		th=0.8, multicores=FALSE, pdf=NULL, maxiterations=10000,
		usebics=FALSE, cores=2, lambda=NULL, B=NULL,Z=Z,maxiter=30) {
	diag(B) <- 0
	# initialise
	#dat[is.na(dat)] <- 0
	antibodies <- rownames(dat)
	tps <- unique(sapply(colnames(dat), function(x) strsplit(x,"_")[[1]][2]))
	reps <- ((ncol(dat)/length(tps))/length(stimuli))
	longprop <- 1:max(length(tps),(nrow(phi)*100))
	gammaposs <- propagate.effect.set(phi,longprop,stimuli,reps=reps)
	gammaposs <- uniquegammaposs(gammaposs)
	# now get an initial gamma matrix
	gammax <- propagate.effect.set(phi,tps,stimuli,reps=reps)
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
	movetypes <- c("switchtype","delete","addactivation","addinhibition","revert")
	bestmodel <- list(phi=phi,L=Linit,aic=aicinit,bic=bicinit,posterior=posteriorinit,dat=dat,
			theta=thetax, gamma=gammax, gammaposs=gammaposs, tps=tps, stimuli=stimuli, reps=reps,
			maxiter=maxiter, TSA=NULL, Tt=NULL, lastmove="addactivation", coords=c(1,1), lambda=lambda,
			B=B,Z=Z,pegm=1,pegmundo=1,nummoves=length(movetypes))
#	phit <- phi
#	for(k in 1:25) {
#		cat(".")
#		posteriorold <- bestmodel$posterior
#		phi <- matrix(sample(c(0,1,2),length(phit),replace=T),nrow=nrow(phit),ncol=ncol(phit),dimnames=dimnames(phit))
#		longprop <- 1:max(length(tps),(nrow(phi)*100))
#		gammaposs <- propagate.effect.set(phi,longprop,stimuli,reps=reps)
#		gammaposs <- uniquegammaposs(gammaposs)
#		# now get an initial gamma matrix
#		gammax <- propagate.effect.set(phi,tps,stimuli,reps=reps)
#		#gammax <- propagate.effect.set(phimax,tps,stimuli,reps=reps)
#		Ltmp <- likl(dat,gammax)
#		Linit <- Ltmp$L
#		thetax <- Ltmp$theta
#		Linit[Linit==Inf] <- 0
#		Linit[Linit==-Inf] <- 0	
#		Linit <- sum(Linit)
#		aicinit <- get.aic(phi,Linit)
#		bicinit <- get.bic(phi,Linit, length(dat))
#		if(is.null(lambda)) {
#			posteriorinit <- NULL
#		} else {
#			posteriorinit <- posterior(phi, sum(Linit), lambda, B, Z)
#		}
#		if(posteriorinit > posteriorold) {
#			print(posteriorinit)
#			bestmodel <- list(phi=phi,L=Linit,aic=aicinit,bic=bicinit,posterior=posteriorinit,dat=dat,
#					theta=thetax, gamma=gammax, gammaposs=gammaposs, tps=tps, stimuli=stimuli, reps=reps,
#					maxiter=maxiter, TSA=NULL, Tt=NULL, lastmove="addactivation", coords=c(1,1), lambda=lambda,
#					B=B,Z=Z)
#		}
#	}	

	it <- 1
	stats <- matrix(0, nrow=maxiterations, ncol=11, dimnames=list(1:maxiterations, c("MAP", "tp","tn","fp","fn","sn","sp","lambda","acpt","lacpt","stmove")))
	freqa <- freqi <- bestmodel$phi
	freqa[freqa!=0] <- 0
	freqi[freqi!=0] <- 0
	while(it < maxiterations) {
		cat("++ ", it, " ")
		newlambda <- runif(1, bestmodel$lambda-1, bestmodel$lambda+1)
		newlambda <- min(max(0.01,newlambda),30)
		movetype <- sample(1:length(movetypes),1)
		if(all(bestmodel$phi==0))
			movetype <- sample(c(3,4),1)
		if(all(bestmodel$phi!=0))
			movetype <- sample(c(2,5),1)
			#movetype <- sample(c(1,2,5),1)
		
		
		st <- system.time(b1 <- mcmc_move(bestmodel, movetypes[movetype]))
		ret <- mcmc_accept(bestmodel, b1, newlambda)
		bestmodel <- ret$bestproposal 
		if(!is.null(phiorig)) {
			O <- phiorig
			M <- bestmodel$phi
			#comp <- compare.graphs.tc(O=phiorig,M=bestmodel$phi)
			#stats[it,] <- as.matrix(cbind(bestmodel$posterior, comp[1:6], bestmodel$lambda))
			tp <- length(which(O==M & (O==1 | O==2)))
			tn <- length(which(O==M & O==0))
			fn <- length(which(O!=M & (O==1 | O==2)))
			fp <- length(which(O!=M & O==0))
			sn <- tp/(tp+fn)
			sp <- tn/(tn+fp)
			stats[it,] <- cbind(bestmodel$posterior, tp, tn, fp, fn, sn, sp, bestmodel$lambda, ret$acpt, ret$lacpt, st[3])
		}
		#if(it>1000) {
			tmp <- bestmodel$phi
			tmp[bestmodel$phi==2] <- 0
			freqa <- freqa + tmp
			tmp <- bestmodel$phi
			tmp[bestmodel$phi==1] <- 0
			freqi <- freqi + (tmp/2)
		#}

		# some convergence statistic

		#SSW <- sd(stats[,"MAP"],na.rm=T)
		#SSB <- 
		#Rhat <- ((nrow(stats)-1)/nrow(stats)) * SSW
		#l <- maxiterations
		l <- it #nrow(stats)
		if(it%%25==1){
			if(!is.null(pdf))
				pdf(pdf)
			par(mfrow=c(4,2),mar=c(3,4,1,1),oma=c(1,1,1,1))
			plot(1:l, stats[1:l,"MAP"], type='l',ylab="",xlab="iteration",main="log posterior")
			plot(1:l, stats[1:l,"sn"], type='l',ylim=c(0,1),ylab="",xlab="iteration",col="blue",main="SN/SP plot")
			lines(1:l, stats[1:l,"sp"], lty=2,col="violet")
			legend("topleft", c("sn","sp"), lty=c(1,2),col=c("blue","violet")) #,"red","orange"))
			plot(1:l, stats[1:l,"acpt"], lty=1, col="#222222",main="acpt")
			plot(1:l, stats[1:l,"lacpt"], lty=1, col="#222222",main="lacpt")
			plot(1:l, stats[1:l,"lambda"], type='l',ylim=c(0,30),ylab="lambda",xlab="iteration")			
			boxplot(as.data.frame(stats[,"stmove"]),ylab="time")
			if(it>1000) {
				boxplot(as.data.frame(stats[(1000:it),c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
						main=paste("avgSN: ", signif(median(stats[(1000:it),"sn"]),digits=4), "avgSP: ", signif(median(stats[(1000:it),"sp"]),digits=4)))
			} else {
				boxplot(as.data.frame(stats[1:it,c("sn","sp","acpt","lacpt")]), ylim=c(0,1),
						main=paste("BURNIN avgSN: ", signif(median(stats[,"sn"]),digits=4), "avgSP: ", signif(median(stats[,"sp"]),digits=4)))
			}
			# partial autocorrelation function:
			R <- acf(stats[1:l,"MAP"])
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
	#bestmodel[["MAPtrace"]] <- maps
	#colnames(stats) <- c("MAP", "tp","tn","fp","fn","sp","sn")
	bestmodel[["stats"]] <- stats
	bestmodel[["freqa"]] <- freqa
	bestmodel[["freqi"]] <- freqi
	bestmodel[["datx"]] <- dat
	bestmodel
}