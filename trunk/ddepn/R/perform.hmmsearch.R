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
### only R
### dat: N x (T x R) matrix, N: number of proteins, T number of timepoints, R: number of replicates
### gammaprime: N x (T x R) matrix
### gammaposs: N x M matrix, M: number of reachable states derived in effect propagation
### E: M x (T x R) matrix: Emission matrix
### A: M x M matrix: Transition matrix
### viterbi: M x T matrix: which path to take
perform.hmmsearch <- function(phi.n, bestmodel) {
	cat(".")
	tps <- bestmodel$tps
	stimuli <- bestmodel$stimuli
	dat <- bestmodel$dat
	hmmiterations <- bestmodel$hmmiterations
	gamprimetotal <- NULL
	gamposstotal <- NULL
	# separate HMM for each experiment, i.e. each stimulus
	for(s in stimuli) {
		exind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(dat))
		R <- length(exind)/length(tps)
		datx <- dat[,exind]
		gammaposs <- propagate.effect.simple(phi.n,stimulus=s,stimuli=stimuli)
		colnames(gammaposs) <- paste(paste(names(s),collapse="&"), colnames(gammaposs), sep="_")
		V <- rownames(datx)
		TC <- unique(colnames(datx))
		M <- ncol(gammaposs)
		N <- nrow(datx)
		T <- length(tps)
		## initial transition matrix
		Adimn <- colnames(gammaposs)
		# all transitions equally likely
		#A <- matrix((1/(M*M)),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		## random transition matrix
		A <- matrix(runif(M*M,0,1),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		A <- A/sum(A,na.rm=TRUE)
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
			Liktmp <- likl(datx,gamprime)
			Lik <- Liktmp$L
			thetaprime <- Liktmp$theta
			Lik[Lik == Inf] <- 0
			Lik[Lik == -Inf] <- 0
			Lik <- sum(Lik,na.rm=TRUE)
			if(Lold!=-Inf) {
				if(abs((abs(Lik)-abs(Lold))) <= 1) 
					break
				diff <- abs((abs(Lik)-abs(Lold)))
				diffsold[it] <- diff
				## check for repeating patterns
				if(diff %in% diffsold[-length(diffsold)]) {
					if(length(which(diffsold==diff))>25) {
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
						# only up to 5 restarts
						if(restarts < 5){
							gamprime <- replicatecolumns(gammaposs[,sort(sample(M,T,replace=TRUE))],R)
							it <- 0
							Lik <- -Inf
							diffold <- -100
							equally <- 0
							next
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
			A <- A - log2(sum(2^A, na.rm=TRUE))	
		} # end while
		# now we have an A, an E and a gammaprime for the first experiment
		# save the gammaprime
		gamprimetotal <- cbind(gamprimetotal, gamprime)
		gamposstotal <- cbind(gamposstotal, gammaposs)
	} # end outer for loop
	Liktmp <- likl(dat,gamprimetotal)
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
			gammaposs=gamposstotal)
	return(L.res)
}
