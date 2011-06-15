# Use either netga or mcmc output lists from ddepn,
# generate a final theta-matrix
# plot the data profiles of each protein and indicate
# the active and passive distributions in the plots.
# 
# Author: benderc
###############################################################################

## use a return list for multiple mcmc runs to get the 
## consensus theta matrix
get_theta_consensus <- function(ret) {
	rets <- ret$samplings
	nr <- nrow(rets[[1]]$theta)
	nc <- ncol(rets[[1]]$theta)
	reps <- length(rets)
	mufinal <- array(NA, dim=c(nr,nc,reps))
	for(i in 1:length(rets)) {
		ret <- rets[[i]]
		## mu_run holds the running mean of the parameters
		## sd_run holds the running sd of the parameters -> don't use it yet
		## think about how to use these deviations
		mutmp <- ret$theta
		#mutmp <- ret$mu_run
		mufinal[,,i] <- mutmp
	}
	thetafin <- apply(mufinal, c(1,2), median, na.rm=TRUE)
	colnames(thetafin) <- colnames(mutmp)
	rownames(thetafin) <- rownames(mutmp)
	thetafin
}
plot_profiles <- function(ret, log=FALSE, ord=NULL,
		mfrow=c(2,2), plothist=FALSE, selection.criterion="aic") {
	if(class(ret)=="list") {
		P <- ret$P
		## mcmc or netga?
		if(is.null(P)) { #mcmc
			dat <- ret$samplings[[1]]$dat
			theta <- get_theta_consensus(ret)
		} else { # netga
			dat <- ret$dat
			thetas <- NULL
			for(i in 1:length(P)) {
				thetas <- cbind(thetas, P[[i]]$theta)
			}
			Parr <- array(thetas, dim=c(nrow(thetas),(ncol(thetas)/length(P)),length(P)))
			theta <- apply(Parr,c(1,2),mean,na.rm=TRUE)
			thetasd <- apply(Parr,c(1,2),sd,na.rm=TRUE)
		}
	} else { # only a data matrix should be plot
		dat <- ret
		theta <- NULL
	}
	stimuli <- which(apply(dat, 1, function(x) all(x==0)))
	if(log) {
		dat <- log2(dat)
	}
	## find all experiments
	expers.fac <- unique(sub("_[0-9]*$","",colnames(dat)))
	## sort the experiments to make the plots easier to compare
	if(!is.null(ord)) {
		expers.fac <- expers.fac[ord]
	}
	expers <- as.character(sapply(colnames(dat), function(x) strsplit(x, "_")[[1]][1]))
	tps <- as.numeric(sapply(colnames(dat), function(x) strsplit(x, "_")[[1]][2]))
	time <- unique(tps)
	reps <- table(expers)/length(time)
	xn <- seq(0,max(time),0.1)
	constants <- NULL
	ddmats <- NULL
	cnt <- 1
	## Spline fits and plots through the data for each protein
	for(j in 1:nrow(dat)) {
		print(rownames(dat)[j])
		if(j %in% stimuli)
			next
		if(cnt %% (mfrow[1]*mfrow[2]) == 1)
			par(mfrow=mfrow) #min(6,length(expers.fac))))
		cnt <- cnt + 1
		ddmat <- NULL
		## Spline fits for each experiment
		for(i in 1:length(expers.fac)) {
			expf <- expers.fac[i]
			cat(expf," ")
			ind <- which(expers==expf)
			y <- dat[j,ind]
			tp <- as.numeric(sapply(colnames(dat)[ind], function(x) strsplit(x,"_")[[1]][2]))
			## sort by time
			ord <- order(tp)
			y <- y[ord]
			tp <- tp[ord]
			## make a fit through the data points 
			if(all(y==-Inf) || all(y==0) || all(is.na(y))){
				pred <- rep(median(dat,na.rm=TRUE), length(xn))
			} else {
				pred <- bestgam(y,tp,xn,selection.criterion=selection.criterion)
				if(all(pred==pred[1])) {
					constants <- c(constants, j)
					#pred <- jitter(pred) # + rnorm(length(pred), pred[1], 0.0001)
				}
			}
			ddmat <- rbind(ddmat, pred)
		}		
		yl <- ifelse(log,"log2","")
		## plot the data for each experiment, teh fits and parameters
		for(i in 1:length(expers.fac)) {
			expf <- expers.fac[i]
			## the splines
			plot(xn, ddmat[i,],type="l",
					ylab=paste(yl,"Intensity"),xlab="time [min]",ylim=range(ddmat),
					main=c(paste("Stimulus:",expf),paste("Protein:",rownames(dat)[j])))
			ind <- which(expers==expf)
			y <- dat[j,ind]
			tp <- as.numeric(sapply(colnames(dat)[ind], function(x) strsplit(x,"_")[[1]][2]))
			mt <- match(time, unique(tp))
			mt <- mt[!is.na(mt)]
			at <- time[mt]
			## the data
			boxplot(y~tp,add=TRUE,at=at,border="#08080850", axes=FALSE)
			## the parameters
			if(!is.null(theta)) {
				abline(h=theta[j,1],col="red")
				abline(h=theta[j,1]+theta[j,2],lty=3,col="red")
				abline(h=theta[j,1]-theta[j,2],lty=3,col="red")
				text(1,theta[j,1],"mu active")
				abline(h=theta[j,3],col="green")
				abline(h=theta[j,3]+theta[j,4],lty=3,col="green")
				abline(h=theta[j,3]-theta[j,4],lty=3,col="green")
				text(1,theta[j,3],"mu passive")
			}
		}
	}
	cnt <- 1
	## should the histograms be plotted?
	if(plothist) {
		for(j in 1:nrow(dat)) {
			if(j %in% stimuli)
				next
			#if(cnt %% 9 == 1)
			#	par(mfrow=mfrow)
			cnt <- cnt + 1
			h <- hist(dat[j,], main=rownames(dat)[j], breaks=20, freq=FALSE)
			if(!is.null(theta)) {
				xx <- seq(min(h$breaks),max(h$breaks),length.out=1000)
				yya <- dnorm(xx,theta[j,1],theta[j,2])
				lines(xx, yya, col="red")
				text(theta[j,1],max(yya),"active",col="red")
				yyp <- dnorm(xx,theta[j,3],theta[j,4])
				lines(xx, yyp,col="green")
				text(theta[j,3],max(yyp),"passive",col="green")
			}
		}
	}
}
