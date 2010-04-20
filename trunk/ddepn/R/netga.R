# genetic algorithm for network search
# n: initial population size
# q: selection rate
# m: mutation rate
# nodes: names of the nodes in the network


netga <- function(datx, stimuli, P=NULL, maxiterations=1000, p=100,
		q=0.3, m=0.8, maxiter=30, multicores=FALSE, usebics=FALSE, cores=2,
		lambda=NULL, B=NULL,
		Z=NULL, scorefile="score.pdf") {
  datx[is.na(datx)] <- 0
  V <- rownames(datx)
  tps <- unique(sapply(colnames(datx), function(x) strsplit(x,"_")[[1]][2]))
  reps <- ((ncol(datx)/length(tps))/length(stimuli))
  phireference <- matrix(0,nrow=length(V), ncol=length(V), dimnames=list(V,V))

  #####################################################################
  #  First create a population of networks if none is given           #
  #####################################################################
  if(is.null(P)) {
	  phireference <- matrix(0,nrow=length(V), ncol=length(V), dimnames=list(V,V))
	  X <- vector("list",p)
	  X[[1]] <- phireference
	  notstim <- setdiff(1:nrow(datx),unlist(stimuli))
	  phireference[,notstim] <- phireference[,notstim] + 1
	  diag(phireference) <- 0
	  X[[2]] <- phireference
	  if(multicores) {
		P <- mclapply(X, getfirstphi, datx=datx, stimuli=stimuli, V=V, tps=tps, reps=reps, maxiter=maxiter, lambda=lambda, B=B, Z=Z, mc.preschedule=FALSE,mc.cores=cores)		
	  } else {
		P <- lapply(X, getfirstphi, datx=datx, stimuli=stimuli, V=V, tps=tps, reps=reps, maxiter=maxiter, lambda=lambda, B=B, Z=Z)
	  }
  }
  if(any(sapply(P, class)!="list")) {
	  print("netga.R: Some elements in the network list P seem to be empty.")
	  browser()
  }
  ## check if all individuals are set correctly
  ## is this a problem of mclapply/lapply? sometimes, returnvalues in the list P are null...
  for(i in 1:length(P)) {
	  if(class(P[[i]])=="try-error" || is.null(P[[i]])){
		  P[[i]] <- getfirstphi(X[[i]], datx=datx, stimuli=stimuli, V=V, tps=tps, reps=reps, maxiter=maxiter, lambda=lambda, B=B, Z=Z)
	  }
  }  
  if(usebics){
	  wks <- sapply(P, function(x) x$bic)
	  optwks <- quantile(wks)[3]
	  oldoptwks <- -Inf
	  wks2 <- wks - max(wks) -1
	  probs <- wks2/sum(wks2)
	  #optwks <- -Inf
  } else {
	  if(is.null(lambda)) {
		  wks <- sapply(P, function(x) x$L)
	  } else {
		  wks <- sapply(P, function(x) x$posterior)
	  }
	  optwks <- quantile(wks)[3]
	  oldoptwks <- Inf
	  wks2 <- wks + abs(min(wks)) +1
	  probs <- wks2/sum(wks2)
	  #optwks <- Inf
  }
  
  #####################
  # main loop         #
  #####################
  p <- length(P)
  diffpercent <- opts <- NULL
  numequalscore <- 0
  autoc <- list(acf=rep(0,5),n.used=10)
  for(i in 1:maxiterations) {
#   #### doesn't seem to be sufficient for convergence, complaint by a reviewer from ismb
#	## test if the mean of the last 20 differences of the average fitness values are 
#	## equal to 0: assume it to be normally distributed, therefore perform a t-test
#	## if H0 cannot be rejected anymore, then we reached convergence
#	if(i>20) {
#		pdiff <- t.test(diff(opts[(length(opts)-20):length(opts)]))$p.value	  
#	} else {
#		pdiff <- 0
#	}
#	if(pdiff>0.1) {
#		for(ii in 1:length(P)) {
#			P[[ii]][["iterations"]] <- i
#		}
#		break
#	}
	pdiff <- "x"
	# terminate criterion: autocorrelation ???
	# if all autocorrelation values for all lags are bigger than 7*sqrt(n)/n, then 
	# stop calculation
	if(all(autoc$acf > 7*sqrt(autoc$n.used)/autoc$n.used)) {
		print("YEEEEEEHAAAA, autocor found stop criterion.")
		for(ii in 1:length(P)) {
			P[[ii]][["iterations"]] <- i
		}
		break
	}
	# terminate criterion: 10x equal optimal score, then return
	if(oldoptwks==optwks) {
		numequalscore <- numequalscore + 1
		print(paste("No improvement in optimal score for ", numequalscore, " times."))
		if(numequalscore==10) {
			for(ii in 1:length(P)) {
				P[[ii]][["iterations"]] <- i
			}
			break
		}
	} else {
		numequalscore <- 0
	}
	##new generation
    Pprime <- list()
	#browser()
    ################
    # selection    #
    ################
	#browser()
	diffpercent <- c(diffpercent,abs(round(((min(wks) - mean(wks))/min(wks)*100), digits=3)))
	opts <- c(opts, optwks)
    print(paste("selection ",i, "diff(opt,avg): ", diffpercent[length(diffpercent)], " Diffs in Opts==0? ", pdiff))
	if(i %% 10 == 0) {
		if(!is.null(pdf)) {
			pdf(scorefile,width=8,height=10)			
		}
		#browser()
		par(mfrow=c(5,1))
		plot(opts, type='l', ylab="median scores", xlab="generation", main=paste("min: ", round(min(opts),3), " max: ", round(max(opts),3)))
		autoc <- acf(opts, main="Autocorrelation of median scores",ci.type="ma")
		abline(h=-2*sqrt(autoc$n.used)/autoc$n.used,col="orange",lty=4)
		abline(h=2*sqrt(autoc$n.used)/autoc$n.used,col="orange",lty=4)
		pautoc <- pacf(opts, main="Partial autocorrelation of median scores")
		abline(h=-2*sqrt(pautoc$n.used)/pautoc$n.used,col="orange",lty=4)
		abline(h=2*sqrt(pautoc$n.used)/pautoc$n.used,col="orange",lty=4)
		plot(diff(opts), type='b', ylab="differences avg scores (t-1 -> t)",xlab="generation i+1", pch="*")
		plot(diffpercent, type='l', ylab="percent optimumscore - avgscore", xlab="generation",main=paste("optimum(minimum difference): ",min(diffpercent)))
		if(!is.null(pdf)) {
			dev.off()	
		}
	}
	#write.table(cbind(diffpercent, opts), file="diffsopts.txt")
    # select only the better models, i.e. maximise the likelihoods
	# or minimize bics
	if(usebics) {
		selection <- which(wks < optwks) # minimise the bic
		if(length(selection)==0)
			selection <- sample(which(wks==min(wks)),1)
	} else {
		selection <- which(wks > optwks) # maximise the L or posterior
		if(length(selection)==0)
			selection <- sample(which(wks==max(wks)),1)
	}
    # if selected too much, then sample from the possibilities
	# with probability proportional to the score
    if(length(selection)>ceiling((p*(1-q)))) {
	  probs2 <- probs[selection]
      selection <- sample(selection, ceiling((p*(1-q))), prob=probs2)
    }
	Pprime <- c(Pprime, P[selection])
    # define barrier that must be exceeded in the next run
	# optimum before crossover and mutation
	if(usebics) {
		print(paste("Selected: ", length(Pprime), " models with bic < ", oldoptwks, ". New minbic: ", optwks))
	} else {
		if(is.null(lambda)) {
			print(paste("Selected: ", length(Pprime), " models with L > ", oldoptwks, ". New maxL: ", optwks))
		} else {
			print(paste("Selected: ", length(Pprime), " models with post > ", oldoptwks, ". New maxPosterior: ", optwks))
		}
	}
	oldoptwks <- optwks
	# get the number of crossovers
    # size of population - already selected individuals
    # i.e. perform crossover for all not selected individuals
    numcrossings <- p - length(selection)
    numcrossings <- (numcrossings - numcrossings%%2)

    ##############
    # crossover  #
    ##############
    print(paste("crossover",i))
	# sample randomly number of crossings individuals, i.e. numcrossings/2 pairs
	# take preferrably the most fit individuals for crossingover
    crossing <- matrix(sample(1:p, numcrossings, prob=probs, replace=FALSE), nrow=2)
	# sampling with uniform prob
	#crossing <- matrix(sample(1:p, numcrossings, replace=FALSE), nrow=2)
	phicrosses <- list()
	PP <- list()
	# do all crossovers
	counter <- 1
    for(k in 1:ncol(crossing)) {
      phicross <- crossover(P[[crossing[1,k]]]$phi, P[[crossing[2,k]]]$phi)
	  PP[[counter]] <- P[[crossing[1,k]]]
	  PP[[counter]]$phi <- phicross[[1]]
	  PP[[(counter+1)]] <- P[[crossing[2,k]]]
	  PP[[(counter+1)]]$phi <- phicross[[2]]
	  counter <- counter+2
	}
#browser()
#	if(any(sapply(PP, class)!="list")) {
#		browser()
#	}

	if(multicores) {
		ret <- mclapply(PP, function(x) perform.hmmsearch(x$phi, x), mc.preschedule=FALSE,mc.cores=cores)	
  	} else {
		ret <- lapply(PP, function(x) perform.hmmsearch(x$phi, x))	
	}
	for(k in 1:length(ret)) {
		if(is.null(ret[[k]])){
			cat("*!!*")
			x <- PP[[k]]
			xret <- perform.hmmsearch(x$phi,x)
			ret[[k]] <- xret
		}
		bestmodel <- PP[[k]]
		L.res <- ret[[k]]
		bestmodel$gamma <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
		bestmodel$theta <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
		bestmodel$L <- L.res$Likl
		bestmodel$bic <- L.res$bic
		bestmodel$aic <- L.res$aic
		if(is.null(lambda)) {
			bestmodel$posterior <- NULL
		} else {
			bestmodel$posterior <- posterior(bestmodel$phi, bestmodel$L, lambda, B, Z)
		}
		bestmodel$gammaposs <- L.res$gammaposs
		Pprime <- c(Pprime, list(bestmodel))
	}

	#Pprime <- c(Pprime, PP)
 	# if new generation is smaller than old one, add some unchanged individuals
    if(length(Pprime)<p) {
      rest <- sample(setdiff(1:p, union(selection, as.vector(crossing))),(p-length(Pprime)))
      Pprime <- c(Pprime, P[rest])
    }
	## now i have changed individuals, get the scores again
	if(usebics) {
		wks <- sapply(Pprime, function(x) x$bic)
		wks2 <- wks - max(wks) -1
		probs <- wks2/sum(wks2)
		#optwks <- quantile(wks,na.rm=T)[2] # 25% quantile
		optwks <- quantile(wks,na.rm=T)[3] # median
	} else {
		# maximize over posterior if prior is given
		if(is.null(lambda)) {
			wks <- sapply(Pprime, function(x) x$L)
		} else {
			wks <- sapply(Pprime, function(x) x$posterior)
		}
		wks2 <- wks + abs(min(wks)) +1
		probs <- wks2/sum(wks2)
		#optwks <- quantile(wks,na.rm=T)[4] # 75% quantile
		optwks <- quantile(wks,na.rm=T)[3] # median
	}
	
    ##############
    # mutation   #
    ##############
	# mutate p*m individuals
    # select the individuals proportional to their likelihood?
    #mutation <- sample(1:p, (p*m), prob=(1-probs)) ## mutate the worst ones
	mutation <- sample(1:p, (p*m), prob=probs) ## mutate the best ones
	#mutation <- sample(1:p, (p*m), prob=(1-probs)) ## mutate with uniform probability
    print(paste("mutation",i))
	counter <- 1
	oldphis <- list()
	oldedges <- newedges <- NULL
	for(k in mutation) {
		phitmp <- Pprime[[k]]$phi
		oldphis[[counter]] <- phitmp
		diag(phitmp) <- -1
		phitmp[,unlist(stimuli)] <- matrix(-1,nrow=nrow(datx),ncol=length(unlist(stimuli)))
		position <- sample(which(phitmp!=-1),1) # select the position
		types <- setdiff(c(0,1,2),phitmp[position])
		type <- sample(types, 1) # select the type of the edge to which it will mutate
		phi.n <- Pprime[[k]]$phi
		oldedges <- c(oldedges, phi.n[position])
		newedges <- c(newedges, type)
		phi.n[position] <- type
		Pprime[[k]]$phi <- phi.n
		counter <- counter + 1
	}
	cat("performing ", length(mutation), " mutations..")
	if(multicores) {
		ret <- mclapply(Pprime[mutation], function(x){perform.hmmsearch(x$phi, x)}, mc.preschedule=FALSE,mc.cores=cores)	
	} else {
			ret <- lapply(Pprime[mutation], function(x) perform.hmmsearch(x$phi, x))	
	}
#	## fang ab, dass manchmal leere Werte zurückgegeben werden, warum???
	for(k in 1:length(ret)) {
		if(is.null(ret[[k]])){
			cat("*!*")
			x <- Pprime[mutation][[k]]
			xret <- perform.hmmsearch(x$phi,x)
			ret[[k]] <- xret
		}
		## here the mutated network is already contained
		bestmodel <- Pprime[[mutation[k]]]
		L.res <- ret[[k]]
		posteriornew <- NULL
		if(usebics) {
			scorenew <- -L.res$bic
			scoreold <- -optwks
		} else {
			if(is.null(lambda)) {
				scorenew <- L.res$Likl
				scoreold <- optwks
			} else {
				posteriornew <- posterior(bestmodel$phi, L.res$Likl, lambda, B, Z)
				scorenew <- posteriornew
				scoreold <- optwks
			}
		}
		# if score is better, than accept the new parameters, state matrix etc.
		if(scorenew > scoreold) {
			print(paste("Improved score by mutation: old edge: ", oldedges[k], " new edge: ",newedges[k], "oldscore: ", scoreold, "newscore: ", scorenew))
			bestmodel$gamma <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
			bestmodel$theta <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
			bestmodel$L <- L.res$Likl
			bestmodel$bic <- L.res$bic
			bestmodel$aic <- L.res$aic
			bestmodel$posterior <- posteriornew
			bestmodel$gammaposs <- L.res$gammaposs
		} else {
			# else restore the old network; parameters and matrices are still the old ones
			bestmodel$phi <- oldphis[[k]]
		}
		Pprime[[mutation[k]]] <- bestmodel
	}
	if(usebics) {
		wks <- sapply(Pprime, function(x) x$bic)
		wks2 <- wks - max(wks) -1
		probs <- wks2/sum(wks2)
		optwks <- quantile(wks,na.rm=T)[3] # median
	} else {
		# maximize over posterior if prior is given
		if(is.null(lambda)) {
			wks <- sapply(Pprime, function(x) x$L)
		} else {
			wks <- sapply(Pprime, function(x) x$posterior)
		}
		wks2 <- wks + abs(min(wks)) +1
		probs <- wks2/sum(wks2)
		optwks <- quantile(wks,na.rm=T)[3] # median
	}
    P <- Pprime
	garbage <- gc(verbose=FALSE)
  } # end main loop
  return(P)
}
