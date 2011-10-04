# genetic algorithm for network search
# p: initial population size
# q: selection rate
# m: mutation rate

netga <- function(dat, stimuli, P=NULL, maxiterations=1000, p=100,
		q=0.3, m=0.8, hmmiterations=30, multicores=FALSE, usebics=FALSE, cores=2,
		lambda=NULL, B=NULL,
		Z=NULL, scorefile=NULL,fanin=4,
		gam=NULL,it=NULL,K=NULL,quantL=.5,quantBIC=.5, priortype="none", plotresults=TRUE,
		scale_lik=FALSE, allow.stim.off=TRUE,debug=0,retobj=NULL) {
  dat[is.na(dat)] <- 0
  V <- rownames(dat)
  ordstim <- sapply(stimuli, function(x) paste(names(x),collapse="&"))
  tmp <- tapply(colnames(dat), gsub("_.*$","",colnames(dat)), get_reps_tps)
  tmp <- tmp[ordstim]	
  tps <- lapply(tmp, function(x) x$tps)
  reps <- sapply(tmp, function(x) as.numeric(x$reps))
  phireference <- matrix(0,nrow=length(V), ncol=length(V), dimnames=list(V,V))

  #####################################################################
  #  First create a population of networks if none is given           #
  #####################################################################
  if(is.null(P)) {
	  X <- vector("list",p)
	  X[[1]] <- phireference
	  if(fanin>=nrow(phireference)) {
	  	notstim <- setdiff(1:nrow(dat),unlist(stimuli))
	  	phireference[,notstim] <- phireference[,notstim] + 1 # set all zeros to 1, i.e. make a completely connected network with only activations
	  	diag(phireference) <- 0 ## no self activations
	  	X[[2]] <- phireference
	  }
	  if(multicores) {
		P <- mclapply(X, getfirstphi, dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,
				                      fanin=fanin,gam=gam,it=it,K=K,priortype=priortype,scale_lik=scale_lik, allow.stim.off=allow.stim.off,
									  mc.preschedule=FALSE,mc.cores=cores)		
	  } else {
		P <- lapply(X, getfirstphi, dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,
				                    fanin=fanin,gam=gam,it=it,K=K,priortype=priortype,scale_lik=scale_lik, allow.stim.off=allow.stim.off)
	  }
  }
  ## check if all individuals are set correctly
  ## is this a problem of mclapply/lapply? sometimes, returnvalues in the list P are null...
  for(i in 1:length(P)) {
	  if(length(which(colSums(ifelse(P[[i]]$phi==0,0,1))>fanin))>0 |
			  any(P[[i]]$phi[,unique(names(unlist(stimuli)))]!=0)) {
		  print("**********  INIT: fanin too big somewhere or edge leading to stimulus")
		  browser()
	  }	  
	  if(class(P[[i]])=="try-error" || is.null(P[[i]])){
		  P[[i]] <- getfirstphi(X[[i]], dat=dat,stimuli=stimuli,V=V,tps=tps,reps=reps,hmmiterations=hmmiterations,lambda=lambda,B=B,Z=Z,fanin=fanin,gam=gam,it=it,K=K,priortype=priortype,scale_lik=scale_lik, allow.stim.off=allow.stim.off)
	  }
	  # store the GA parameters
	  #P[[i]]$q <- q
	  #P[[i]]$m <- m
	  #P[[i]]$usebics <- usebics
	  #P[[i]]$quantBIC <- quantBIC
	  #P[[i]]$quantL <- quantL
  }  
  if(any(sapply(P, class)!="list")) {
	  print("netga.R: Some elements in the network list P seem to be empty.")
	  browser()
  }
  if(usebics){
	  wks <- sapply(P, function(x) x$bic)
	  score_quantile <- quantile(wks,na.rm=T,probs=quantBIC)
	  old_score_quantile <- Inf
	  wks2 <- wks - max(wks) -1
	  probs <- wks2/sum(wks2)
  } else {
	  if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {
		  wks <- sapply(P, function(x) x$posterior)
	  } else {
		  wks <- sapply(P, function(x) x$L)
	  }
	  score_quantile <- quantile(wks,na.rm=T,probs=quantL)
	  old_score_quantile <- -Inf
	  wks2 <- wks + abs(min(wks)) +1
	  probs <- wks2/sum(wks2)
  }
  
  #####################
  # main loop         #
  #####################
  p <- length(P)
  diffpercent <- NULL
  #diffpercent <- opts <- NULL
  numequalscore <- 0
  autoc <- list(acf=rep(0,5),n.used=10)
  start <- 1
  ## define a matrix holding some statistics on the development of the scores
  scorestats <- matrix(NA, nrow=maxiterations, ncol=18)
  colnames(scorestats) <- c("dL_total","dP_total",
			  "dL_crossover","dL_mutation",
			  "dP_crossover","dP_mutation",
			  "dL_total_abs","dP_total_abs",
			  "dL_crossover_abs","dL_mutation_abs",
			  "dP_crossover_abs","dP_mutation_abs",
			  "liklihood","prior",
			  "liklihood_mad","prior_mad",
			  "score", "score_mad")
  rownames(scorestats) <- 1:maxiterations
  ## if retobj is given, the GA should be resumed
  if(!is.null(retobj)) {
	start <- nrow(retobj$scorestats)+1
	scorestats[1:(start-1),] <- retobj$scorestats
  }
  #####################
  ### GA Main loop
  #####################
  cls()
  if(debug==0) {
	print("###########################")
	print("# Using genetic algorithm #")
	print("###########################")
  	pb <- txtProgressBar(min = 0, max = maxiterations, style = 3)
  }
  for(iter in start:maxiterations) {
	if(debug==0) {
		#cls() # clear the screen
		#pcnt <- round(iter/maxiterations * 100)
		#print(paste("[",pcnt, "%] done.",sep=""))
		setTxtProgressBar(pb, iter)
		
	}
	#pdiff <- "x"
	# terminate criterion: 50x equal optimal score, then return
	if(old_score_quantile==score_quantile) {
		numequalscore <- numequalscore + 1
		print(paste("No improvement in optimal score for ", numequalscore, " times."))
		if(numequalscore==50) {
			for(ii in 1:length(P)) {
				P[[ii]][["iterations"]] <- iter
			}
			break
		}
	} else {
		numequalscore <- 0
	}
	##new generation
    Pprime <- list()
    ################
    # selection    #
    ################
	diffpercent <- c(diffpercent,abs(round(((min(wks) - mean(wks))/min(wks)*100), digits=3)))
	#opts <- c(opts, score_quantile)
	if((debug==1 & iter%%50==1) | debug==2) { # print info depending on debug level
		print(paste("Iteration ",iter))
    	print(paste("   Selection step. Difference optimum to average: ", diffpercent[length(diffpercent)]))
	}
	#########################
	### plot some population diagnostics
	if(iter %% 10 == 1 & iter>1) {
		if(plotresults) {	
			if(!is.null(scorefile)) {
				pdf(scorefile,width=8,height=10)			
			}
			layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
			## score trace
			opts <- scorestats[1:(iter-1),"score"]
			plot(opts, type='l', ylab="median scores", xlab="generation", main=paste("Score trace: min: ", round(min(opts),2), " max: ", round(max(opts),2)))
			plot(diff(opts), type='l', ylab="difference median scores", xlab="generation", main=paste("Median score diffs; min:",round(min(diff(opts)),digits=2),"max:",round(min(diff(opts)),digits=2)))
			## liklihood trace
			opts <- scorestats[1:(iter-1),"liklihood"]
			plot(opts, type='l', ylab="liklihood", xlab="generation", main=paste("Liklihood trace: min: ", round(min(opts),2), " max: ", round(max(opts),2)))
			plot(diff(opts), type='l', ylab="difference liklihood", xlab="generation", main=paste("Median liklihood diffs; min:",round(min(diff(opts)),digits=2),"max:",round(min(diff(opts)),digits=2)))
			## prior trace
			opts <- scorestats[1:(iter-1),"prior"]
			plot(opts, type='l', ylab="prior", xlab="generation", main=paste("Prior trace: min: ", round(min(opts),2), " max: ", round(max(opts),2)))
			plot(diff(opts), type='l', ylab="difference liklihood", xlab="generation", main=paste("Median prior diffs; min:",round(min(diff(opts)),digits=2),"max:",round(min(diff(opts)),digits=2)))	
	#		layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE))
	#		plot(opts, type='l', ylab="median scores", xlab="generation", main=paste("Median score trace: min: ", round(min(opts),3), " max: ", round(max(opts),3)))
	#		autoc <- acf(opts, main="Autocorrelation of median scores",ci.type="ma")
	#		abline(h=-2*sqrt(autoc$n.used)/autoc$n.used,col="orange",lty=4)
	#		abline(h=2*sqrt(autoc$n.used)/autoc$n.used,col="orange",lty=4)
	#		pautoc <- pacf(opts, main="Partial autocorrelation of median scores")
	#		abline(h=-2*sqrt(pautoc$n.used)/pautoc$n.used,col="orange",lty=4)
	#		abline(h=2*sqrt(pautoc$n.used)/pautoc$n.used,col="orange",lty=4)
	#		plot(diff(opts), type='b', ylab="differences avg scores (t-1 -> t)",xlab="generation i+1", pch="*")
	#		plot(diffpercent, type='l', ylab="percent optimumscore - avgscore", xlab="generation",main=paste("optimum(minimum difference): ",min(diffpercent)))
			if(!is.null(scorefile)) {
				dev.off()	
			}
		}
	}
    # select only the better models, i.e. maximise the likelihoods
	# or minimize bics
	if(usebics) {
		selection <- which(wks < score_quantile) # minimise the bic
		# don't stop if no score is less than the quantile of the scores
		if(length(selection)==0)
			selection <- sample(which(wks==min(wks)),1)
	} else {
		selection <- which(wks > score_quantile) # maximise the L or posterior
		# don't stop if no score is less than the quantile of the scores
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
	if((debug==1 & iter%%50==1) | debug==2) { # print info depending on debug level
		if(usebics) {
			print(paste("      Selected: ", length(Pprime), " models with bic < ", old_score_quantile, ". New minbic: ", score_quantile))
		} else {
			if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {
			#if(laplace || scalefree) {
				print(paste("      Selected: ", length(Pprime), " models with post > ", old_score_quantile, ". New maxPosterior: ", score_quantile))
			} else {
				print(paste("      Selected: ", length(Pprime), " models with L > ", old_score_quantile, ". New maxL: ", score_quantile))
			}
		}
	}
	old_score_quantile <- score_quantile
	# get the number of crossovers
    # size of population - already selected individuals
    numcrossings <- p - length(selection)
    numcrossings <- (numcrossings - numcrossings%%2)
	################
    ################
    ## crossover  ##
    ################
	################
	if((debug==1 & iter%%50==1) | debug==2) { # print info depending on debug level
    	print(paste("   Crossover step."))
	}
	# sample randomly number of crossings individuals, i.e. numcrossings/2 pairs
	# take preferrably the most fit individuals for crossingover
    crossing <- matrix(sample(1:p, numcrossings, prob=probs, replace=FALSE), nrow=2)
	# sampling with uniform prob
	phicrosses <- list()
	PP <- list()
	# do all crossovers
	counter <- 1
    for(k in 1:ncol(crossing)) {
      phicross <- crossover(P[[crossing[1,k]]]$phi, P[[crossing[2,k]]]$phi,fanin=fanin)
	  if(length(which(colSums(ifelse(phicross[[1]]==0,0,1))>fanin))>0 |
			  any(phicross[[1]][,unique(names(unlist(stimuli)))]!=0)) {
		  print("**********  CROSS1: fanin too big somewhere or edge leading to stimulus")
		  browser()
	  }
	  if(length(which(colSums(ifelse(phicross[[2]]==0,0,1))>fanin))>0 ||
	  			any(phicross[[2]][,unique(names(unlist(stimuli)))]!=0)) {
		  print("**********   CROSS2: fanin too big somewhere or edge leading to stimulus")
		  browser()
	  }
  
	  PP[[counter]] <- P[[crossing[1,k]]]
	  PP[[counter]]$phi <- phicross[[1]]
	  PP[[(counter+1)]] <- P[[crossing[2,k]]]
	  PP[[(counter+1)]]$phi <- phicross[[2]]
	  counter <- counter+2
	}
	if(multicores) {
		ret <- mclapply(PP, function(x) perform.hmmsearch(x$phi, x), mc.preschedule=FALSE,mc.cores=cores)	
  	} else {
		ret <- lapply(PP, function(x) perform.hmmsearch(x$phi, x))	
	}
	## init statistics vector
	dL <- dP <- NULL
	for(k in 1:length(ret)) {
		if(is.null(ret[[k]])){
			cat("*!!*")
			x <- PP[[k]]
			xret <- perform.hmmsearch(x$phi,x)
			ret[[k]] <- xret
		}
		bestmodel <- PP[[k]] ## contains already the new network
		L.res <- ret[[k]]
		## statistics, Liklihood difference 
		dL <- c(dL,bestmodel$L - L.res$Likl)
		bestmodel$gamma <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
		bestmodel$theta <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
		bestmodel$L <- L.res$Likl
		bestmodel$bic <- L.res$bic
		bestmodel$aic <- L.res$aic
		if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {
			## new prior; L.res$phix should be equal to bestmodel$phi
			if(!all(L.res$phix==bestmodel$phi)) {
				print("Error: New network proposal is not set correctly.")
				browser()
			}
			prnew <- prior(bestmodel$phi, lambda, B, Z, gam, it, K, priortype)
			## statistics, prior difference
			dP <- c(dP,bestmodel$pr - prnew)
			## setting the new prior and posterior
			bestmodel$pr <- prnew
			bestmodel$posterior <- bestmodel$L + bestmodel$pr
		} else {
			## statistics, prior difference
			dP <- c(dP,0)
			prnew <- 0
			bestmodel$pr <- 0
			bestmodel$posterior <- 0
		}
		bestmodel$gammaposs <- L.res$gammaposs
		Pprime <- c(Pprime, list(bestmodel))
	}
	## bestmodel: unchanged model, L.res holds the scores for the changed model
	scorestats[iter,"dL_crossover"] <- median(dL)
	scorestats[iter,"dL_crossover_abs"] <- median(abs(dL))
	## statistics: old prior - new prior
	scorestats[iter,"dP_crossover"] <- median(dP)
	scorestats[iter,"dP_crossover_abs"] <- median(abs(dP))
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
		score_quantile <- quantile(wks,na.rm=T,probs=quantBIC)
	} else {
		# maximize over posterior if prior is given
		if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {	
			wks <- sapply(Pprime, function(x) x$posterior)
		} else {
			wks <- sapply(Pprime, function(x) x$L)
		}
		wks2 <- wks + abs(min(wks)) +1
		probs <- wks2/sum(wks2)
		score_quantile <- quantile(wks,na.rm=T,probs=quantL)
	}
	
    ##############
    # mutation   #
    ##############
	# mutate p*m individuals
    # select the individuals proportional to their likelihood?
    #mutation <- sample(1:p, (p*m), prob=(1-probs)) ## mutate the worst ones
	mutation <- sample(1:p, (p*m), prob=probs) ## mutate the best ones
	if((debug==1 & iter%%50==1) | debug==2) { # print info depending on debug level
    	print(paste("   Mutation step."))
	}
	counter <- 1
	oldphis <- list()
	oldedges <- newedges <- NULL
	for(k in mutation) {
		phitmp <- Pprime[[k]]$phi
		oldphis[[counter]] <- phitmp
		# take care of fanin
		out <- which(colSums(ifelse(phitmp==0,0,1))>=fanin)
		diag(phitmp) <- -1	
		phitmp[,unique(names(unlist(stimuli)))] <- matrix(-1,nrow=nrow(dat),ncol=length(unique(unlist(stimuli))))
		fout <- which(which(phitmp==0,arr.ind=TRUE)[,2]%in%out)
		phitmp[which(phitmp==0)[fout]] <- -1
		position <- sample(which(phitmp!=-1),1) # select the position
		types <- setdiff(c(0,1,2),phitmp[position])
		type <- sample(types, 1) # select the type of the edge to which it will mutate
		phi.n <- Pprime[[k]]$phi
		oldedges <- c(oldedges, phi.n[position])
		newedges <- c(newedges, type)
		phi.n[position] <- type
		if(length(which(colSums(ifelse(phi.n==0,0,1))>fanin))>0 ||
				any(phi.n[,unique(names(unlist(stimuli)))]!=0)) {
			print("********** MUTATION! fanin too big somewhere or edge leading to stimulus...")
			browser()
		}
		Pprime[[k]]$phi <- phi.n
		counter <- counter + 1
	}
	if(debug==2) { # print info depending on debug level
		print(paste("      Performing ", length(mutation), " mutations."))
	}
	if(multicores) {
		ret <- mclapply(Pprime[mutation], function(x){perform.hmmsearch(x$phi, x)}, mc.preschedule=FALSE,mc.cores=cores)	
	} else {
		ret <- lapply(Pprime[mutation], function(x) perform.hmmsearch(x$phi, x))	
	}
	## init statistics vector
	dLm <- dPm <- NULL
	## fang ab, dass manchmal leere Werte zurueckgegeben werden, warum???
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
		## statistics, Likelihood difference
		dLm <- c(dLm,bestmodel$L - L.res$Likl)
		if(usebics) {
			scorenew <- -L.res$bic
			scoreold <- -score_quantile
			dPm <- c(dPm,0)
			prnew <- 0
			posteriornew <- 0
		} else {
			if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform") && !usebics) {
				prnew <- prior(bestmodel$phi, lambda, B, Z, gam, it, K, priortype )
				## statistics: prior difference
				dPm <- c(dPm,bestmodel$pr - prnew)
				posteriornew <- L.res$Likl + prnew
				scorenew <- posteriornew
				scoreold <- score_quantile
			} else {
				## statistics: prior difference
				dPm <- c(dPm,0)
				prnew <- 0
				posteriornew <- 0
				scorenew <- L.res$Likl
				scoreold <- score_quantile
			}
		}
		# if score is better, than accept the new parameters, state matrix etc.
		if(scorenew > scoreold) {
			#if((debug==1 & iter%%50==1) | debug==2) { # print info depending on debug level
			if(debug==2) { # print info depending on debug level
				#print(paste("Improved score by mutation: old edge: ", oldedges[k], " new edge: ",newedges[k], "oldscore: ", scoreold, "newscore: ", scorenew))
				print(paste("   Improved score: oldscore: ", scoreold, "newscore: ", scorenew))
			}
			bestmodel$gamma <- matrix(L.res$gammax,nrow=nrow(bestmodel$gamma),ncol=ncol(bestmodel$gamma),dimnames=dimnames(bestmodel$gamma))
			bestmodel$theta <- matrix(L.res$thetax,nrow=nrow(bestmodel$theta),ncol=ncol(bestmodel$theta),dimnames=dimnames(bestmodel$theta))
			bestmodel$L <- L.res$Likl
			bestmodel$bic <- L.res$bic
			bestmodel$aic <- L.res$aic
			bestmodel$posterior <- posteriornew
			bestmodel$pr <- prnew
			bestmodel$gammaposs <- L.res$gammaposs
		} else {
			# else restore the old network; parameters and matrices are still the old ones
			bestmodel$phi <- oldphis[[k]]
		}
		Pprime[[mutation[k]]] <- bestmodel
	}
	## bestmodel: unchanged model, L.res holds the scores for the changed model
	scorestats[iter,"dL_mutation"] <- median(dLm)
	scorestats[iter,"dL_mutation_abs"] <- median(abs(dLm))
	scorestats[iter,"dL_total"] <- median(c(dLm,dL))
	scorestats[iter,"dL_total_abs"] <- median(abs(c(dLm,dL)))
	## statistics: old prior - new prior
	scorestats[iter,"dP_mutation"] <- median(dPm)
	scorestats[iter,"dP_mutation_abs"] <- median(abs(dPm))
	scorestats[iter,"dP_total"] <- median(c(dPm,dP))
	scorestats[iter,"dP_total_abs"] <- median(abs(c(dPm,dP)))
	## statistics: L, Pr, score 
	scorestats[iter,"liklihood"] <- median(sapply(Pprime, function(x) x$L))
	scorestats[iter,"prior"] <- median(sapply(Pprime, function(x) x$pr))
	scorestats[iter,"liklihood_mad"] <- mad(sapply(Pprime, function(x) x$L))
	scorestats[iter,"prior_mad"] <- mad(sapply(Pprime, function(x) x$pr))

	if(usebics) {
		wks <- sapply(Pprime, function(x) x$bic)
		wks2 <- wks - max(wks) -1
		probs <- wks2/sum(wks2)
		score_quantile <- quantile(wks,na.rm=T,probs=quantBIC)
	} else {
		# maximize over posterior if prior is given
		if(priortype %in% c("laplaceinhib","laplace","scalefree","uniform")) {		
			wks <- sapply(Pprime, function(x) x$posterior)	
		} else {
			wks <- sapply(Pprime, function(x) x$L)
		}
		wks2 <- wks + abs(min(wks)) +1
		probs <- wks2/sum(wks2)
		score_quantile <- quantile(wks,na.rm=T,probs=quantL)
	}
	## statistics: score
	scorestats[iter,"score"] <- score_quantile
	scorestats[iter,"score_mad"] <- mad(wks)
    P <- Pprime
	garbage <- gc(verbose=FALSE)
	if(debug==1)
		cat(".")
  } # end main loop		
  if(debug==0) {	  
 	 close(pb)
	 print("done.")
  }
  return(list(P=P,scorestats=scorestats))
}
