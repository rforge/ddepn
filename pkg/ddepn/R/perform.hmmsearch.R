# TODO: Add comment
# 
# Author: benderc
###############################################################################
## E-step
## get emission probabilities
getE <- function(x, datx, thetaprime) {
	L <- log2(pmin(dnorm(datx, mean=thetaprime[,"mu.passive"],thetaprime[,"sd.passive"])*ifelse(x==0,1,NA),
					dnorm(datx, mean=thetaprime[,"mu.active"],thetaprime[,"sd.active"])*ifelse(x==1,1,NA),na.rm=TRUE))
	L[L==Inf] <- 0
	L[L==-Inf] <- 0
	L <- colSums(L, na.rm=TRUE)
	L
}

replicatecolumns <- function(mat, replicates=4) {
	coln <- rep(colnames(mat), each=replicates)
	mm <- matrix(apply(mat, 2, rep, times=replicates),nrow=nrow(mat))
	colnames(mm) <- coln
	rownames(mm) <- rownames(mat)
	mm
}

perform.hmmsearch <- function(phi.n, bestmodel) {
	## select implementation
	env <- get("envDDEPN")
    implementation <- get("IMPLEMENTATION", envir=env)
	#implementation <- get("IMPLEMENTATION", pos=globalenv())
	## do this only until different number of timepoints  
	## in each experiment is implemented in C
	#if(class(bestmodel$tps)=="list")
	#	bestmodel$tps <- bestmodel$tps[[1]]
	
	switch(implementation,
			R=perform.hmmsearch_R(phi.n, bestmodel),
			R_globalest=perform.hmmsearch_globalest(phi.n, bestmodel),
			C_globalest=perform.hmmsearch_C_globalest(phi.n, bestmodel),
			ret<-perform.hmmsearch_C(phi.n, bestmodel))
}

### only R
### dat: N x (T x R) matrix, N: number of proteins, T number of timepoints, R: number of replicates
### gammaprime: N x (T x R) matrix
### gammaposs: N x M matrix, M: number of reachable states derived in effect propagation
### E: M x (T x R) matrix: Emission matrix
### A: M x M matrix: Transition matrix
### viterbi: M x T matrix: which path to take
perform.hmmsearch_R <- function(phi.n, bestmodel) {
	#cat(".")
	#browser()
	tpsall <- bestmodel$tps
	tps <- unlist(tpsall) ## vector representation for c-function
	Tall <- sapply(tpsall, length)
	stimuli <- bestmodel$stimuli
	dat <- bestmodel$dat
	hmmiterations <- bestmodel$hmmiterations
	scale_lik <- bestmodel$scale_lik
	allow.stim.off <- bestmodel$allow.stim.off
	gamprimetotal <- NULL
	gamposstotal <- NULL
	# separate HMM for each experiment, i.e. each stimulus
	#for(s in stimuli) {
	for(i in 1:length(stimuli)) {
		s <- stimuli[[i]]
		tps <- tpsall[[i]]
		T <- Tall[i]
		exind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(dat))
		#R <- length(exind)/length(tps)
		R <- bestmodel$reps[i]
		datx <- dat[,exind]
		gammaposs <- propagate.effect.simple(phi.n,stimulus=s,stimuli=stimuli,allow.stim.off=allow.stim.off)
		colnames(gammaposs) <- paste(paste(names(s),collapse="&"), colnames(gammaposs), sep="_")
		V <- rownames(datx)
		TC <- unique(colnames(datx))
		M <- ncol(gammaposs)
		N <- nrow(datx)
		#T <- length(tps)
		## initial transition matrix
		Adimn <- colnames(gammaposs)
		# all transitions equally likely
		#A <- matrix((1/(M*M)),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		## random transition matrix
		A <- matrix(runif(M*M,0,1),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		A <- A/rowSums(A,na.rm=TRUE) ## valid transition matrix has row sums of 1
		#A <- A/sum(A,na.rm=TRUE)
		pseudocount <- 1
		pseudocountsum <- M
		A <- log2(A)
		## initial gamprime 
		gamprime <- replicatecolumns(gammaposs[,sort(sample(M,T,replace=TRUE))],R)
		## initial theta
		Lik <- -Inf
		diffold <- -100
		diffsold <- rep(NA, hmmiterations)
		equally <- 0
		it <- 0
		restarts <- 0		
		while(it <= hmmiterations) {
			it <- it + 1
			## total likelihood
			Lold <- sum(Lik,na.rm=TRUE)
			Liktmp <- likl(datx,gamprime,scale_lik=scale_lik)
			Lik <- Liktmp$L
			thetaprime <- Liktmp$theta
			Lik[Lik == Inf] <- 0
			Lik[Lik == -Inf] <- 0
			Lik <- sum(Lik,na.rm=TRUE)
			#print(paste("Lik:",Lik,"Lold:",Lold,"Diff Lik Lold: ",abs((abs(Lik)-abs(Lold)))))
			if(Lold!=-Inf) {
				if(abs((abs(Lik)-abs(Lold))) <= 1) 
					break
				diff <- abs((abs(Lik)-abs(Lold)))
				diffsold[it] <- diff
				## check for repeating patterns
				if(diff %in% diffsold[-length(diffsold)]) {
					if(length(which(diffsold==diff))>10) {
						equally <- 2
						diff <- diffold
					}
				}
				if(diff==diffold) {
					equally <- equally + 1
					# restart if switching behaviour occurs,
					# don't know where this comes from
					if(equally == 3) {
						restarts <- restarts + 1
						# only up to a number of restarts
						if(restarts < 3){
							gamprime <- replicatecolumns(gammaposs[,sort(sample(M,T,replace=TRUE))],R)
							#it <- 0
							Lik <- -Inf
							diffold <- -100
							equally <- 0
							diffsold <- rep(NA, hmmiterations)
							next
						} else {
							break
						}
					}		
				}
				diffold <- diff
			}
			E <- t(apply(gammaposs, 2, getE, datx=datx, thetaprime=thetaprime))
			dimnames(E) <- list(colnames(gammaposs),colnames(datx))
			# transform to 3d array
			sel <- rep(seq(1,T*R,by=R),R)
			seladd <- rep(seq(0,R-1),each=T)
			sel <- sel + seladd
			E <- array(E[,sel,drop=F],dim=c(M,T,R))		
			## viterbi maximum likelihood path		
			viterbi <- matrix(0,nrow=M,ncol=T,dimnames=list(colnames(gammaposs),unique(colnames(datx))))
			maxima <- rep(-Inf,T)
			maxima.ind <- rep(0,T)
			al <- -log2(M)
			# j is position in datamatrix
			for(j in 1:ncol(E)) {
				# i ist state to which we jump
				## inner loop
				if(j == 1) {
					vtcol <- rowSums(al + E[,j,,drop=FALSE]) 
				} else {
					vtcol <- apply(A, 2, function(Acol, viterbi, j) max(viterbi[,j-1,drop=FALSE] + Acol, na.rm=TRUE), viterbi=viterbi, j=j) + rowSums(E[,j,,drop=FALSE])
				}
				viterbi[,j] <- vtcol
				vt <- max(vtcol)
				i <- which(vtcol==vt)[1]
				if(maxima[j]<vt) {
					maxima[j] <- vt
					maxima.ind[j] <- i # this is the optimal state series
				}
			}
			maxima <- rep(maxima, each=R)
			maxima.ind <- rep(maxima.ind, each=R)
			## get new gamma suggestion and estimate parameters
			gamprime <- gammaposs[,maxima.ind]
			colnames(gamprime) <- colnames(datx)
			## M-step
			## maximize the transition matrices parameters
			sel <- cbind(1:length(maxima.ind),2:(length(maxima.ind)+1))
			## transitions hold the switchings from state i to state i+1 in the state sequence
			transitions <- matrix(maxima.ind[sel],ncol=2)
			## change the probabilities for transitions that were observerd
			trans <- table(transitions[-nrow(transitions),1])
			transall <- table(paste(transitions[-nrow(transitions),1], transitions[-nrow(transitions),2], sep="_"))		
			ind <- match(as.numeric(sapply(names(transall), function(x) strsplit(x, split="_")[[1]][1])),as.numeric(names(trans)))	
			transprob <- log2((transall+pseudocount)) - log2(trans[ind]+pseudocountsum)
			indices <- sapply(names(transprob), function(x,rows) (as.numeric(strsplit(x, "_")[[1]])-c(0,1)) %*% c(1,rows),rows=nrow(A))
			A[indices] <- transprob
			A <- A - log2(rowSums(2^A, na.rm=TRUE))
			#A <- A - log2(sum(2^A, na.rm=TRUE))	
		} # end while
		# now we have an A, an E and a gammaprime for the first experiment
		# save the gammaprime
		gamprimetotal <- cbind(gamprimetotal, gamprime)
		gamposstotal <- cbind(gamposstotal, gammaposs)
	} # end outer for loop
	#print(paste("ITERATIONS PERFORMED:", it, " / ", hmmiterations))
	Liktmp <- likl(dat,gamprimetotal,scale_lik=scale_lik)
	Lik <- Liktmp$L
	thetaprime <- Liktmp$theta
	Lik[Lik == Inf] <- 0
	Lik[Lik == -Inf] <- 0
	Lik <- sum(Lik,na.rm=TRUE)
	aic <- get.aic(phi.n, Lik)
	bic <- get.bic(phi.n, Lik, length(dat))
	L.res <- list(datx=dat, phix=phi.n, stimx=stimuli,
			gammax=gamprimetotal, thetax=thetaprime,
			replicates=R, Likl=Lik,aic=aic,bic=bic,
			statespace_maxiterations=hmmiterations, 
			gammaposs=gamposstotal,scale_lik=scale_lik,allow.stim.off=allow.stim.off)
	return(L.res)
}

#### uses c-library
#### dat: N x (T x R) matrix, N: number of proteins, T number of timepoints, R: number of replicates
#### gammaprime: N x (T x R) matrix
#### gammaposs: N x M matrix, M: number of reachable states derived in effect propagation
#### E: M x (T x R) matrix: Emission matrix
#### A: M x M matrix: Transition matrix
#### viterbi: M x T matrix: which path to take
perform.hmmsearch_C <- function(phi.n, bestmodel) {
	#browser()
	#cat(".")
	tpsall <- bestmodel$tps
	tps <- unlist(tpsall) ## vector representation for c-function
	T <- sapply(tpsall, length)
	#tps <- bestmodel$tps
	#T <- length(tps)
	stimuli <- bestmodel$stimuli
	dat <- bestmodel$dat
	hmmiterations <- bestmodel$hmmiterations
	GS <- matrix(0, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
	gammaposs <- propagate.effect.set(phi.n, stimuli=stimuli)
	G <- gammaposs
	#G[G!=0] <- 0
	TH <- matrix(0.0, nrow=nrow(dat), ncol=4, dimnames=dimnames(bestmodel$theta))
	stimids <- unlist(stimuli)
	stimnames <- names(stimids)
	stimgrps <- sapply(stimuli, length)
	numexperiments <- length(bestmodel$reps)
	## get the M_i: the number of possible states for each experiment, same dimension of R
	Ms <- table(gsub("_[0-9]$","",colnames(gammaposs)))
	## sort according to the experiments, as they are in the data matrices
	Ms <- Ms[order(match(names(Ms), sub("_[0-9]$","",colnames(gammaposs))))]
	Lik <- 0.0

	# separate HMM for each experiment, i.e. each stimulus	
	ret <- .C("perform_hmmsearch",P=as.integer(phi.n), N=as.integer(nrow(dat)),
			T=as.integer(T), R=as.integer(bestmodel$reps), X=as.double(dat),
			GS=as.integer(GS), G=as.integer(G), Glen=as.integer(length(G)),
			TH=as.double(TH), 	tps=as.integer(tps), stimids=as.integer(stimids-1),
			stimgrps=as.integer(stimgrps), numexperiments=as.integer(numexperiments),
			Likx=as.double(Lik), hmmiterations=as.integer(bestmodel$hmmiterations),
			Msx=as.integer(Ms), PACKAGE="ddepn", NAOK=TRUE)
	
	Lik <- ret$Likx
	aic <- get.aic(phi.n, Lik)
	bic <- get.bic(phi.n, Lik, length(dat))
	gamprimetotal <- matrix(ret$GS, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
	thetaprime <- matrix(ret$TH, nrow=nrow(TH), ncol=ncol(TH), dimnames=dimnames(TH))
	gamposstotal <- matrix(ret$G, nrow=nrow(dat), ncol=ncol(G), dimnames=dimnames(G))
	
	L.res <- list(datx=dat, phix=phi.n, stimx=stimuli,
			gammax=gamprimetotal, thetax=thetaprime,
			replicates=bestmodel$reps, Likl=Lik,aic=aic,bic=bic,
			statespace_maxiterations=hmmiterations, 
			gammaposs=gamposstotal)
	return(L.res)
}

### only R
### !!!estimate theta using all experiments together, do not assemble the total gammprime matrix
### and theta matrix afterwards!!!
### dat: N x (T x R) matrix, N: number of proteins, T number of timepoints, R: number of replicates
### gammaprime: N x (T x R) matrix
### gammaposs: N x M matrix, M: number of reachable states derived in effect propagation
### E: M x (T x R) matrix: Emission matrix
### A: M x M matrix: Transition matrix
### viterbi: M x T matrix: which path to take
## TODO: (1) rename use_C argument in ddepn to implementation, with 3 arguments:
## 1) "R_globalest": use the global parameter estimation over all experiments (possibly slower than the other implementations)
## 2) "R": use the R-version with separate parameter estimation for each experiment
## 3) "C": use the C-version with separate parameter estimation for each experiment
## TODO: (2) add argument "C_globalest": C version of 1)
perform.hmmsearch_globalest <- function(phi.n, bestmodel) {
	#browser()
	#cat(".")
	tpsall <- bestmodel$tps
	#T <- length(tps)
	stimuli <- bestmodel$stimuli
	dat <- bestmodel$dat
	cols <- colnames(dat)
	#reps <- table(sub("_[0-9].*$","",cols)) / length(tps)
	reps <- bestmodel$reps
	hmmiterations <- bestmodel$hmmiterations
	scale_lik <- bestmodel$scale_lik
	allow.stim.off <- bestmodel$allow.stim.off
	gamprimetotal <- NULL
	gamposstotal <- NULL
	# separate HMM for each experiment, i.e. each stimulus
	## first, get the state matrices for each experiment
	initgamp <- NULL
	Mold <- NULL
	Ms <- 0
	for(si in 1:length(stimuli)) {
		s <- stimuli[[si]]
		R <- reps[[si]]
		tps <- tpsall[[si]]
		T <- length(tps)
		exind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(dat))
		#R <- length(exind)/length(tps)
		datx <- dat[,exind]
		gammaposs <- propagate.effect.simple(phi.n,stimulus=s,stimuli=stimuli,allow.stim.off=allow.stim.off)
		colnames(gammaposs) <- paste(paste(names(s),collapse="&"), colnames(gammaposs), sep="_")
		gamposstotal <- cbind(gamposstotal, gammaposs)
		M <- ncol(gammaposs)
		initgamp <- c(initgamp, rep(sort(sample(M,T,replace=TRUE)+Ms),each=R))
		Mold <- c(Mold, M)
		Ms <- Ms + M
	}
	
	V <- rownames(datx)
	TC <- unique(colnames(datx))
	M <- ncol(gamposstotal)
	N <- nrow(datx)
	## initial transition matrix
	Adimn <- colnames(gamposstotal)
	# all transitions equally likely
	#A <- matrix((1/(M*M)),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
	## random transition matrix
	A <- matrix(runif(M*M,0,1),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
	## only the submatrices on the diagonal should have values, 
	## transition from e.g. EGF_1 to EGF&HRG_1 should not be possible
	for(s in stimuli) {
		sel <- which(gsub("_[0-9]$","",rownames(A))==paste(names(s), collapse="&"))
		A[sel,-sel] <- 0
		A[-sel,sel] <- 0
	}
	A <- A/rowSums(A,na.rm=TRUE) ## valid transition matrix has row sums of 1
	#A <- A/sum(A,na.rm=TRUE)
	pseudocount <- 1
	pseudocountsum <- M
	A <- log2(A)
	## initial gamprime
	gamprime <- gamposstotal[,initgamp]
	## initial theta
	Lik <- -Inf
	diffold <- -100
	diffsold <- rep(NA, hmmiterations)
	equally <- 0
	it <- 0
	restarts <- 0	
	#gammaposs <- gamposstotal
	while(it <= hmmiterations) {
		it <- it + 1
		## total likelihood
		Lold <- sum(Lik,na.rm=TRUE)
		Liktmp <- likl(dat,gamprime,scale_lik=scale_lik)
		Lik <- Liktmp$L
		thetaprime <- Liktmp$theta
		Lik[Lik == Inf] <- 0
		Lik[Lik == -Inf] <- 0
		Lik <- sum(Lik,na.rm=TRUE)
		#print(paste("Lik:",Lik,"Lold:",Lold,"Diff Lik Lold: ",abs((abs(Lik)-abs(Lold)))))
		if(Lold!=-Inf) {
			if(abs((abs(Lik)-abs(Lold))) <= 1) 
				break
			diff <- abs((abs(Lik)-abs(Lold)))
			diffsold[it] <- diff
			## check for repeating patterns
			if(diff %in% diffsold[-length(diffsold)]) {
				if(length(which(diffsold==diff))>10) {
					equally <- 2
					diff <- diffold
				}
			}
			if(diff==diffold) {
				equally <- equally + 1
				# restart if switching behaviour occurs,
				# don't know where this comes from
				if(equally == 3) {
					restarts <- restarts + 1
					# only up to a number of restarts
					if(restarts < 3){
						initgamp <- NULL
						Ms <- 0									
						for(si in 1:length(Mold)) {
							tps <- tpsall[[si]]
							T <- length(tps)
							M <- Mold[si]
							#print(M)
							R <- reps[si]
							initgamp <- c(initgamp, rep(sort(sample(M,T,replace=TRUE)+Ms),each=R))
							Ms <- Ms + M
						}
						gamprime <- gamposstotal[,initgamp]
						#gamprime <- replicatecolumns(gammaposs[,sort(sample(M,T,replace=TRUE))],R)
						#it <- 0
						Lik <- -Inf
						diffold <- -100
						equally <- 0
						next
					} else {
						break
					}
				}		
			}
			diffold <- diff
		}
		## hmm for each stimulus separately
		gamprime <- NULL
		An <- A
		#for(s in stimuli) {
		for(si in 1:length(stimuli)) {
			s <- stimuli[[si]]
			tps <- tpsall[[si]]
			T <- length(tps)
			R <- reps[[si]]
			toswitch <- which(!apply(thetaprime, 1, function(x) any(is.na(x))|all(x==0)))
			thetaprime.backup <- thetaprime
			done <- FALSE
			## consistency loop
			## if inconsitencies are found, then modify the thetaprime matrix
			## the idea is that inconsitencies are found when active and passive 
			## states are switched (whatever is active should be passive and vice 
			## versa). Since the likelihood will be equal in both cases, we have 
			## to use the structure of the network to conclude on which intensity
			## levels belong to which state
			##
			## runtime will increase a lot when doing it for all possible substitutions
			## in thetaprime, so just try one substitution of a random number of 
			## switchable rows in thetaprime
			exind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(dat))
			gpind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(gamposstotal))
			Aind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(A))
			datx <- dat[,exind,drop=FALSE]
			gammaposs <- gamposstotal[,gpind,drop=FALSE]
			M <- ncol(gammaposs)
			Asel <- A[Aind,Aind,drop=FALSE] 
			while(TRUE) {
				## viterbi algorithm
				maxima.ind <- viterbi(gammaposs, datx, thetaprime, T, R, M, Asel)
				## get new gamma suggestion and estimate parameters
				tmp <- gammaposs[,maxima.ind]
				colnames(tmp) <- colnames(datx)
				ic <- is_consistent(phi.n, tmp)
				## found inconsitencies?
				## switch some of the active/passive state parameter values
				if(ic$inc>0 && length(toswitch)>0 && !done) {
					#print(paste("Inconsitensies in state series. Repeat HMM with modified thetaprime."))
					thtmp <- thetaprime.backup
					nselect <- sample(1:length(toswitch),1) # number of rows to switch
					l2sw <- sample(toswitch,nselect,replace=FALSE)
					thtmp[l2sw,c(1,2)] <- thetaprime.backup[l2sw, c(3,4)]
					thtmp[l2sw,c(3,4)] <- thetaprime.backup[l2sw, c(1,2)]
					thetaprime <- thtmp
					done <- TRUE
				} else {
					break
				}
			}
			gamprime <- cbind(gamprime,tmp)
			Asel <- updateA(maxima.ind, pseudocount, M, Asel)
			An[Aind,Aind] <- Asel
		}
		A <- An
		A <- A - log2(rowSums(2^A, na.rm=TRUE))
	} # end while
	#print(paste("ITERATIONS PERFORMED:", it, " / ", hmmiterations))
	Liktmp <- likl(dat,gamprime,scale_lik=scale_lik)
	Lik <- Liktmp$L
	thetaprime <- Liktmp$theta
	Lik[Lik == Inf] <- 0
	Lik[Lik == -Inf] <- 0
	Lik <- sum(Lik,na.rm=TRUE)
	aic <- get.aic(phi.n, Lik)
	bic <- get.bic(phi.n, Lik, length(dat))
	L.res <- list(datx=dat, phix=phi.n, stimx=stimuli,
			gammax=gamprime, thetax=thetaprime,
			replicates=R, Likl=Lik,aic=aic,bic=bic,
			statespace_maxiterations=hmmiterations, 
			gammaposs=gamposstotal,scale_lik=scale_lik,allow.stim.off=allow.stim.off)
	return(L.res)
}

#### uses c-library with global parameter estimation
#### dat: N x (T x R) matrix, N: number of proteins, T number of timepoints, R: number of replicates
#### gammaprime: N x (T x R) matrix
#### gammaposs: N x M matrix, M: number of reachable states derived in effect propagation
#### E: M x (T x R) matrix: Emission matrix
#### A: M x M matrix: Transition matrix
#### viterbi: M x T matrix: which path to take
perform.hmmsearch_C_globalest <- function(phi.n, bestmodel) {
	#cat(".")
	tpsall <- bestmodel$tps
	tps <- unlist(tpsall) ## vector representation for c-function
	T <- sapply(tpsall, length)	
	#tps <- bestmodel$tps
	#T <- length(tps)
	stimuli <- bestmodel$stimuli
	dat <- bestmodel$dat
	hmmiterations <- bestmodel$hmmiterations
	GS <- matrix(0, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
	gammaposs <- propagate.effect.set(phi.n, stimuli=stimuli)
	G <- gammaposs
	TH <- matrix(0.0, nrow=nrow(dat), ncol=4, dimnames=dimnames(bestmodel$theta))
	stimids <- unlist(stimuli)
	stimnames <- names(stimids)
	stimgrps <- sapply(stimuli, length)
	numexperiments <- length(bestmodel$reps)
	Lik <- 0
	## get the M_i: the number of possible states for each experiment, same dimension of R
	Ms <- table(gsub("_[0-9]$","",colnames(gammaposs)))
	## sort according to the experiments, as they are in the data matrices
	Ms <- Ms[order(match(names(Ms), sub("_[0-9]$","",colnames(gammaposs))))]
		
	# separate HMM for each experiment, i.e. each stimulus	
	ret <- .C("perform_hmmsearch_globalest",P=as.integer(phi.n), N=as.integer(nrow(dat)),
			T=as.integer(T), R=as.integer(bestmodel$reps), X=as.double(dat),
			GS=as.integer(GS), G=as.integer(G), Glen=as.integer(length(G)),
			TH=as.double(TH), 	tps=as.integer(tps), stimids=as.integer(stimids-1),
			stimgrps=as.integer(stimgrps), numexperiments=as.integer(numexperiments),
			Likx=as.double(Lik), hmmiterations=as.integer(hmmiterations), 
			Msx=as.integer(Ms), PACKAGE="ddepn", NAOK=TRUE)
	
	Lik <- ret$Likx
	aic <- get.aic(phi.n, Lik)
	bic <- get.bic(phi.n, Lik, length(dat))
	gamprimetotal <- matrix(ret$GS, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
	thetaprime <- matrix(ret$TH, nrow=nrow(TH), ncol=ncol(TH), dimnames=dimnames(TH))
	gamposstotal <- matrix(ret$G, nrow=nrow(dat), ncol=ncol(G), dimnames=dimnames(G))
	
	L.res <- list(datx=dat, phix=phi.n, stimx=stimuli,
			gammax=gamprimetotal, thetax=thetaprime,
			replicates=bestmodel$reps, Likl=Lik,aic=aic,bic=bic,
			statespace_maxiterations=hmmiterations, 
			gammaposs=gamposstotal)
	return(L.res)
}

## update the transition matrix
updateA <- function(maxima.ind, pseudocount, M, A) {
	sel <- cbind(1:length(maxima.ind),2:(length(maxima.ind)+1))
	## transitions hold the switchings from state i to state i+1 in the state sequence
	transitions <- matrix(maxima.ind[sel],ncol=2)
	## change the probabilities for transitions that were observerd
	trans <- table(transitions[-nrow(transitions),1])
	transall <- table(paste(transitions[-nrow(transitions),1], transitions[-nrow(transitions),2], sep="_"))		
	ind <- match(as.numeric(sapply(names(transall), function(x) strsplit(x, split="_")[[1]][1])),as.numeric(names(trans)))	
	transprob <- log2((transall+pseudocount)) - log2(trans[ind]+M)
	indices <- sapply(names(transprob), function(x,rows) (as.numeric(strsplit(x, "_")[[1]])-c(0,1)) %*% c(1,rows),rows=nrow(A))
	A[indices] <- transprob
	A
}

## the viterbi algorithm
viterbi <- function(gammaposs, datx, thetaprime, T, R, M, A) {
	E <- t(apply(gammaposs, 2, getE, datx=datx, thetaprime=thetaprime))
	dimnames(E) <- list(colnames(gammaposs),colnames(datx))
	# transform to 3d array
	sel <- rep(seq(1,T*R,by=R),R)
	seladd <- rep(seq(0,R-1),each=T)
	sel <- sel + seladd
	E <- array(E[,sel,drop=FALSE],dim=c(M,T,R))		
	viterbi <- matrix(0,nrow=M,ncol=T,dimnames=list(colnames(gammaposs),unique(colnames(datx))))
	maxima <- rep(-Inf,T)
	maxima.ind <- rep(0,T)
	al <- -log2(M)
	# j is position in datamatrix
	for(j in 1:ncol(E)) {
		# i ist state to which we jump
		## inner loop
		if(j == 1) {
			vtcol <- rowSums(al + E[,j,,drop=FALSE]) 
		} else {
			vtcol <- apply(A, 2, function(Acol, viterbi, j) max(viterbi[,j-1,drop=FALSE] + Acol, na.rm=TRUE), viterbi=viterbi, j=j) + rowSums(E[,j,,drop=FALSE])
		}
		viterbi[,j] <- vtcol
		vt <- max(vtcol)
		i <- which(vtcol==vt)[1]
		if(maxima[j]<vt) {
			maxima[j] <- vt
				maxima.ind[j] <- i # this is the optimal state series
		}
	}
	maxima <- rep(maxima, each=R)
	maxima.ind <- rep(maxima.ind, each=R)
	maxima.ind
}

## traverse from leaves up to the root
## find inconsitency markers for each child depending on its parents
##   inconsistency is found whenever:
##   edge pa --> ch
##        gammax[pa,t_i] == 1 # if the parent is active, then state switch from active to passive should not occur at the child
##			gammax[ch, t_i] == 1 & gammax[ch, t_(i+1)] == 0 not possible
##        gammax[pa,t_i] == 0 # if the parent is passive, then state switch from passive to active should not occur at the child
##			gammax[ch, t_i] == 0 & gammax[ch, t_(i+1)] == 1 not possible
##   edge pa --| ch
##        gammax[pa,t_i] == 1 # if parent is active, then state switch from passive to active should not occur
##			gammax[ch, t_i] == 0 & gammax[ch, t_(i+1)] == 1 not possible
##        gammax[pa,t_i] == 0
##			gammax[ch, t_i] == 1 & gammax[ch, t_(i+1)] == 0 not possible
#is_consistent2 <- function(phi, gammax) {
#	phix <- phi
#	inc <- 0
#	allchecks <- allpos <- list()
#	cnt <- 0
#	while(nrow(phix)>1){
#		cnt <- cnt + 1
#		## get the leafs
#		leafs <- which(apply(phix, 1, function(x) all(x==0)))
#		checks <- rep(1, ncol(gammax)-1)
#		checkstab <- NULL
#		for(i in 1:length(leafs)) {
#			leaf <- leafs[i]
#			parents <- which(phix[,leaf]!=0)
#			if(length(parents)>0) {
#				for(j in 1:length(parents)) {
#					pa <- parents[j]
#					consis <- checkC(phix, gammax, leaf, pa)
#					checks <- pmin(checks, consis)
#					checkstab <- cbind(checkstab, consis)
#					#inc <- inc + checkC(phix, gammax, leaf, pa)	
#				}
#			} else {
#				checks <- rep(1, ncol(gammax)-1) ## if no parents, all consistent
#				checkstab <- cbind(checkstab, checks)
#			}
#		}
#		allchecks[[cnt]] <- checkstab
#		allpos[[cnt]] <- which(checks==0)
#		#inc <- inc + sum(ifelse(checks==0,0,1))
#		inc <- inc + sum(1-checks) #ifelse(checks==0,0,1))
#		phix <- phix[-leafs,-leafs,drop=FALSE]
#	}
#	list(inc=inc, allchecks=allchecks, allpos=allpos)
#}
is_consistent <- function(phi, gammax) {
	phix <- phi
	inc <- 0
	allchecks <- allpos <- list()
	for(i in 1:nrow(phix)) {
		nd <- rownames(phix)[i]
		## get the parents
		parents <- which(phix[,nd]!=0)
		checks <- rep(1, ncol(gammax)-1)
		checkstab <- NULL
		if(length(parents)>0) {
			for(j in 1:length(parents)) {
				pa <- parents[j]
				consis <- checkC(phix, gammax, nd, pa)
				checks <- pmin(checks, consis, na.rm=TRUE)
				checks[checks==-14] <- 1 ## special case
				checkstab <- cbind(checkstab, consis)
			}
		} else {
			checks <- rep(1, ncol(gammax)-1) ## if no parents, all consistent
			checkstab <- cbind(checkstab, checks)
		}
		colnames(checkstab) <- rownames(phix)[parents]
		allchecks[[i]] <- checkstab
		allpos[[i]] <- which(checks==0)
		inc <- inc + sum(1-checks) #ifelse(checks==0,0,1))
	}
	names(allchecks) <- names(allpos) <- rownames(phix)
	list(inc=inc, allchecks=allchecks, allpos=allpos)
}

## checks the consistency of a child given one parent and returns the 
## number of inconsistencies
checkC <- function(phi, gammax, ch, pa) {
	ed <- phi[pa,ch]
	nc <- ncol(gammax)
	gpa <- gammax[pa, -nc]
	gch <- gammax[ch, -nc]
	gch_shift <- gammax[ch, 2:nc] # shift left the child state vector
	test <- paste(ed, gpa, gch, gch_shift, sep="_")
	levs <- c("1_0_0_0","1_0_0_1","1_0_1_0","1_0_1_1",
			  "1_1_0_0","1_1_0_1","1_1_1_0","1_1_1_1",
			  "2_0_0_0","2_0_0_1","2_0_1_0","2_0_1_1",
			  "2_1_0_0","2_1_0_1","2_1_1_0","2_1_1_1")
	## consistent represented by value 1
	## inconsistent by val 0, undetermined by NA
	## forced consistencies:
	##    1->0 if pa inhib is 1, marked as special case with -14
	##      this is the only case that can override the inconsistency
	##      at transition 1->0, for active activation parent
	##      
	## forced inconsistencies:
	##    0->1 if pa inhib is 1
	##    1->1 if pa inhib is 1
	vals <- c(1, NA, 1, 1,
			  1, 1, 0, 1,
			  1, 1, 1, 1,
			  1, 0, -14, 0)
	testfac <- factor(test, ordered=TRUE, levels=levs)
	vals[as.numeric(testfac)]		
}


