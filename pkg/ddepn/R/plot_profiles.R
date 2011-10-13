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
		retX <- rets[[i]]
		## mu_run holds the running mean of the parameters
		## sd_run holds the running sd of the parameters -> don't use it yet
		## think about how to use these deviations
		mutmp <- retX$theta
		#mutmp <- ret$mu_run
		mufinal[,,i] <- mutmp
	}
	thetafin <- apply(mufinal, c(1,2), median, na.rm=TRUE)
	sdfin <- apply(mufinal, c(1,2), mad, na.rm=TRUE)
	colnames(thetafin) <- colnames(sdfin) <- colnames(mutmp)
	rownames(thetafin) <- rownames(sdfin) <- rownames(mutmp)
	list(theta=thetafin,sd=sdfin)
}

get_gamma_consensus <- function(ret) {
	rets <- ret$samplings
	nr <- nrow(rets[[1]]$gamma)
	nc <- ncol(rets[[1]]$gamma)
	reps <- length(rets)
	mufinal <- array(NA, dim=c(nr,nc,reps))
	for(i in 1:length(rets)) {
		mufinal[,,i] <- rets[[i]]$gamma
	}
	
	gfin <- apply(mufinal, c(1,2), sum, na.rm=TRUE)
	colnames(gfin) <- colnames(rets[[i]]$dat)
	rownames(gfin) <- rownames(rets[[i]]$dat)
	gfin
}

heatmapcolors <- function(dat,ncol,lowcol="green",highcol="red",middlecol="white") {
	# define color palette, set the 0 to 'black' if possible
	# if all values negative, only take blue values, if all positive, only take red values
	# generates ncol color steps
	rng <- range(dat,na.rm=T)
	if(rng[1]<0){
		# all negative
		if(rng[2]<0) {
			border=ncol
		} else {
			border <- round(abs(rng[1])/abs(rng[2]-rng[1])*ncol)
		}
	} else { # all positive
		border <- 0
	}
	colors <- c(colorpanel(border,low=lowcol,high=middlecol),colorpanel(ncol-border,low=middlecol,high=highcol))
	return(colors)
}

plot_profiles <- function(ret, log=FALSE, ord=NULL,
		mfrow=c(4,4), plotcurves=TRUE, plothist=TRUE, selection.criterion="aic") {
	if(class(ret)=="list") {
		P <- ret$P
		actprof <- NULL
		## mcmc or netga?
		if(is.null(P)) { #mcmc
			dat <- ret$samplings[[1]]$dat
			tmp <- get_theta_consensus(ret)
			theta <- tmp$theta
			thsd <- tmp$sd
			## get a consensus gamma
			gammax <- get_gamma_consensus(ret)
			## get an activity profile from gammax: sum up the numbers of 
			## activities in each row, the higher the number, the more consistent
			## the activity is
			actprof <- apply(gammax, 1, function(xx) tapply(xx, names(xx), sum))
			actprof <- t(actprof[match(unique(colnames(gammax)),rownames(actprof)),])
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
	acol <- "red"
	pcol <- "blue"
	#actprof <- actprof/length(ret$samplings) ## divide by numer of runs
	#actprof <- actprof/ ## divide by number of replicates
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
#browser()
		if(plotcurves) {
			## plot the data for each experiment, the fits and parameters
			for(i in 1:length(expers.fac)) {
				if(!is.null(actprof)) {
					actp <- actprof/reps[i]
					#actcols <- heatmapcolors(actprof, 5, lowcol="green", highcol="red", middlecol="green")
					numcols <- length(ret$samplings)+1 #reps[i]+1
					actcols <- heatmapcolors(actp, numcols, lowcol="green", highcol="red", middlecol="green")
					#actcols <- heatmapcolors((actprof - (max(actprof)/2)), 5, lowcol="green", highcol="red", middlecol="white")
					#actcols <- heatmapcolors(t(scale(t(scale(actprof)))), 5, lowcol="green", highcol="red", middlecol="white")
				} else {
					actcols <- rep("white", reps[i]+1)
				}
				#pad <- rep(actcols[length(actcols)], (reps[i]-5))
				#actcols <- c(actcols, pad) 
				expf <- expers.fac[i]
				actprofexpers <- gsub("_[0-9]*$","",colnames(actp))
				actprofind <- which(actprofexpers==expf)
				## the splines
				plot(xn, ddmat[i,],type="n",
						ylab=paste(yl,"Intensity"),xlab="time [min]",ylim=range(dat[j,], na.rm=TRUE),
						main=c(paste("Stimulus:",expf),paste("Protein:",rownames(dat)[j])))
				ind <- which(expers==expf)
				y <- dat[j,ind]
				tp <- as.numeric(sapply(colnames(dat)[ind], function(x) strsplit(x,"_")[[1]][2]))
				mt <- match(time, unique(tp))
				mt <- mt[!is.na(mt)]
				at <- time[mt]
				## the data
				bcol <- actcols[actp[j,actprofind]+1] #/reps[i])+1)]
				#bcol <- actcols[actprof[j,actprofind]+1]
				boxplot(y~tp,add=TRUE,at=at,border="#08080850", axes=FALSE,col=bcol)
				lines(xn, ddmat[i,],lwd=2)
				## the parameters
				if(!is.null(theta)) {
					abline(h=theta[j,1],col=acol,lwd=2)
					abline(h=theta[j,1]+theta[j,2],lty=3,lwd=1,col=acol)
					abline(h=theta[j,1]-theta[j,2],lty=3,lwd=1,col=acol)
					text(1,theta[j,1],"mu active",col=acol)
					abline(h=theta[j,3],col=pcol,lwd=2)
					abline(h=theta[j,3]+theta[j,4],lty=3,lwd=1,col=pcol)
					abline(h=theta[j,3]-theta[j,4],lty=3,lwd=1,col=pcol)
					text(1,theta[j,3],"mu passive",col=pcol)
				}
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
				lines(xx, yya, col=acol,lwd=1.5)
				text(theta[j,1],max(yya),"active",col=acol)
				yyp <- dnorm(xx,theta[j,3],theta[j,4])
				lines(xx, yyp,col=pcol,lwd=1.5)
				text(theta[j,3],max(yyp),"passive",col=pcol)
			}
		}
	}
}
