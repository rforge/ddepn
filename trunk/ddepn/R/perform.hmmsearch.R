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
	L <- colSums(L, na.rm=T)
	L
}
replicatecolumns <- function(mat, replicates=4) {
	coln <- rep(colnames(mat), each=replicates)
	mm <- matrix(apply(mat, 2, rep, times=replicates),nrow=nrow(mat))
	colnames(mm) <- coln
	rownames(mm) <- rownames(mat)
	mm
}

### datx: N x T x R matrix, N: number of proteins, T number of timepoints, R: number of replicates
### gammaprime: N x T x R matrix
### gammaposs: N x M matrix, M: number of reachable states derived in effect propagation
### E: M x T x R matrix: Emission matrix
### A: M x M matrix: Transition matrix
### viterbi: M x T x R matrix: which path to take
perform.hmmsearch <- function(phi.n, bestmodel) {
	cat(".")
	tps <- bestmodel$tps
	stimuli <- bestmodel$stimuli
	dat <- bestmodel$dat
	maxiter <- bestmodel$maxiter
	R <- bestmodel$reps
	gamprimetotal <- NULL
	gamposstotal <- NULL
	# separates HMM für jedes experiment, also jeden stimulus
	for(s in stimuli) {
		exind <- grep(paste("^",paste(names(s), collapse="&"),"_[0-9]*$",sep=""),colnames(dat))
		datx <- dat[,exind]
		longprop <- 1:max(length(tps),(nrow(phi.n)*100)) # set high maximum number of propagation steps
		gammaposs <- uniquegammaposs(propagate.effect.set(phi.n,longprop,list(s),reps=R))
	#print(paste("*************#################***************"))
	#print(paste("Diff longprop/ncol(gammaposs): ", max(longprop), " / ", ncol(gammaposs)))
	#print(paste("*************#################***************"))
		V <- rownames(datx)
		TC <- unique(colnames(datx))
		M <- ncol(gammaposs)
		N <- nrow(datx)
		T <- length(tps)	
		## initial transition matrix, all transitions equally likely
		Adimn <- colnames(gammaposs) #c(colnames(gamposs),"X0")
		A <- matrix((1/(M*M)),nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		#A <- matrix(0,nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		#A <- apply(A, c(1,2), function(x,A) 1/(ncol(A)^2 - sum(1:(ncol(A)-1))), A=A)
		#A[lower.tri(A)] <- 0
		pseudocount <- 1
		pseudocountsum <- M
		A <- log2(A)
		## pseudocounts brauch ich die??? 
		rkl <- matrix(1,nrow=M,ncol=M,dimnames=list(Adimn,Adimn))
		## initial gamprime: dimensions: 
		#gamprime <- array(rep(gammaposs[,sort(sample(M,T,replace=TRUE))],R), dim=c(N, T, R), dimnames=list(V, TC, 1:R))
		gamprime <- replicatecolumns(gammaposs[,sort(sample(M,T,replace=TRUE))],R)
		## initial theta
		Lik <- -Inf
		diffold <- -100
		equally <- 0
		it <- 0
		restarts <- 0		
		while(it <= maxiter) {
			it <- it + 1
			#print(paste("####################   ",it))
			## total likelihood
			Lold <- sum(Lik,na.rm=T)
			Liktmp <- likl(datx,gamprime)
			Lik <- Liktmp$L
			thetaprime <- Liktmp$theta
			Lik[Lik == Inf] <- 0
			Lik[Lik == -Inf] <- 0
			Lik <- sum(Lik,na.rm=T)
			if(Lold!=-Inf) {
				if(abs((abs(Lik)-abs(Lold))) <= 1) 
					break
				diff <- abs((abs(Lik)-abs(Lold)))
				if(diff==diffold) {
					equally <- equally + 1
					# restart if switching behaviour occurs,
					# don't know where this comes from
					if(equally == 3) {
						#print(paste("####Switching between two local optima, restart.Diff: ",diff))
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
						print("maximum number of restarts reached. Finishing the search.")
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
					vtcol <- apply(A, 2, function(Acol, viterbi, j) max(viterbi[,j-1,drop=FALSE] + Acol, na.rm=T), viterbi=viterbi, j=j) + rowSums(E[,j,,drop=FALSE])
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
			## transitions hold the switchings from state i to state j in the state sequence
			transitions <- matrix(maxima.ind[sel],ncol=2)
			## hier ändere ich nur die übergangswahrscheinlichkeiten, die auch gesehen wurden
			trans <- table(transitions[-nrow(transitions),1])
			transall <- table(paste(transitions[-nrow(transitions),1], transitions[-nrow(transitions),2], sep="_"))		
			ind <- match(as.numeric(sapply(names(transall), function(x) strsplit(x, split="_")[[1]][1])),as.numeric(names(trans)))	
			transprob <- log2((transall+pseudocount)) - log2(trans[ind]+pseudocountsum)
			indices <- sapply(names(transprob), function(x,rows) (as.numeric(strsplit(x, "_")[[1]])-c(0,1)) %*% c(1,rows),rows=nrow(A))
			A[indices] <- transprob
			#A[lower.tri(A)] <- -Inf	
			A <- A - log2(sum(2^A))	
		} # end outer for loop
		# now we have an A, an E and a gammaprime for the first experiment
		# save the gammaprime
		gamprimetotal <- cbind(gamprimetotal, gamprime)
		gamposstotal <- cbind(gamposstotal, gammaposs)
	}
	Liktmp <- likl(dat,gamprimetotal)
	Lik <- Liktmp$L
	thetaprime <- Liktmp$theta
	Lik[Lik == Inf] <- 0
	Lik[Lik == -Inf] <- 0
	Lik <- sum(Lik,na.rm=T)
	aic <- get.aic(phi.n, Lik)
	bic <- get.bic(phi.n, Lik, length(dat))
	L.res <- list(datx=dat, phix=phi.n, stimx=stimuli,
			gammax=gamprimetotal, thetax=thetaprime,
			replicates=R, Likl=Lik,aic=aic,bic=bic,
			statespace_maxiterations=maxiter, 
			gammaposs=gamposstotal)
	return(L.res)
}
