# TODO: Add comment
# 
# Author: benderc
###############################################################################
## get the numbers of the replicates and the timepoints for each experiment
## returns a list containing a list tps with the timepoint labels for each experiment
##         and a vector with the replicate numbers for each experiment
get_reps_tps <- function(nx) {	
	tpsx <- unique(sapply(nx, function(xx) strsplit(xx,"_")[[1]][2]))
	r1 <- table(sub("_[0-9].*$","",nx)) / length(tpsx)
	list(tps=tpsx, reps=r1)
}

## make a pad of NA columns, whenever there are unequal numbers of replicates 
pad_data <- function(dat) {
	tcd <- table(colnames(dat))
	topad <- which(tcd!=max(tcd))
	if(length(topad)>0) {
		xxx <- max(tcd) - tcd[topad]
		cn_pad <- rep(names(xxx), times=xxx)
		pad <- matrix(NA, nrow=nrow(dat), ncol=length(cn_pad), dimnames=list(rownames(dat), cn_pad))
		dat <- cbind(dat, pad)
	}
	dat
}

## sort a data matrix according to experiment, timepoint and replicate indicator
order_experiments <- function(dat) {
	#exper <- gsub("-[0-9]_[0-9]*$","",colnames(dat))
	#brepl <- as.numeric(gsub("(^.*-)([0-9])(_[0-9]*$)","\\2",colnames(dat)))
	#timep <- as.numeric(gsub("^.*_","",colnames(dat)))
	#ord <- order(exper, timep, brepl)
	s1 <- gsub("_.*$","",colnames(dat))
	s2 <- gsub("^.*_","",colnames(dat))
	ord <- order(s1, as.numeric(s2))
	dat <- dat[,ord]
	dat
}


reverse.direction <- function(phi, i, switchtype=FALSE) {
	cs <- coord(i,phi)
	phi[cs[2],cs[1]] <- phi[cs[1],cs[2]]
	phi[cs[1],cs[2]] <- 0
	if(switchtype)
		phi[cs[2],cs[1]] <- ifelse(phi[cs[2],cs[1]]==2,1,2)
	phi
}

## called from simulate.data or simulate.gamma.data.pairs
get.data <- function(gammax,mu.bg=0, sd.bg=0.1,
		mu.signal.a=1, sd.signal.a=0.5,
		mu.signal.i=-1, sd.signal.i=0.5,
		stimulus=NULL,TT=10,R.t=4,R.b=3) {
	
	reps <- R.t*R.b
	datx <- matrix(NA,nrow=nrow(gammax),ncol=ncol(gammax),dimnames=dimnames(gammax))
	nullstate <- matrix(rnorm(nrow(datx)*reps,mu.bg,sd.bg),nrow=nrow(datx),ncol=reps,dimnames=list(rownames(datx),rep(0,reps)))
	nullstate.gamma <- matrix(rep(0,nrow(datx)*reps),nrow=nrow(datx),ncol=reps,dimnames=list(rownames(datx),rep(0,reps)))
	
	datx <- cbind(nullstate,datx)
	gammax <- cbind(nullstate.gamma, gammax)
	tseq <- seq(1,(reps*(TT+1)),by=reps)[-1]
	ss <- 1:nrow(gammax)
	downreg.n <- sample(ss,1) 
	downreg <- sample(ss,downreg.n) # which nodes downreg
	for(t in tseq) {
		for(p in 1:nrow(gammax)) {
			if(p %in% downreg) {
				fac <- -1
			} else {
				fac <- 1
			}
			# stays inactive	
			if(gammax[p,t-1]==0 && gammax[p,t]==0) {
				#datx[p,t:(t+reps-1)] <- rnorm(reps,mu.bg,sd.bg) * fac
				datx[p,t:(t+reps-1)] <- rnorm(reps,mu.bg,sd.bg)
			}
			# becomes active
			if(gammax[p,t-1]==0 && gammax[p,t]==1) {
				#datx[p,t:(t+reps-1)] <- datx[p,(t-reps):(t-1)] + rnorm(reps,mu.signal.a,sd.signal.a) * fac
				datx[p,t:(t+reps-1)] <- rnorm(reps,mu.signal.a,sd.signal.a)
			}
			# becomes inactive, i.e. was inhibited
			if(gammax[p,t-1]==1 && gammax[p,t]==0) {
				#datx[p,t:(t+reps-1)] <- datx[p,(t-reps):(t-1)] - rnorm(reps,mu.signal.i,sd.signal.i) * fac
				#datx[p,t:(t+reps-1)] <- rnorm(reps,mu.signal.i,sd.signal.i)
				datx[p,t:(t+reps-1)] <- rnorm(reps,mu.bg,sd.bg)
			}
			# stays active
			if(gammax[p,t-1]==1 && gammax[p,t]==1) {
				#datx[p,t:(t+reps-1)] <- datx[p,(t-reps):(t-1)] + rnorm(reps,(mu.signal.a/t),sd.signal.a) * fac
				#datx[p,t:(t+reps-1)] <- datx[p,(t-reps):(t-1)] + rnorm(reps,mu.bg,sd.bg)
				datx[p,t:(t+reps-1)] <- rnorm(reps,mu.signal.a,sd.signal.a)
			}
		}
	}
	return(list(datx=datx[,-c(1:reps)],downreg=downreg))
}


# write a matrix to a plot region
plotmatrix <- function(mat,name="") {
	plot.new() # defines new plot with x/y region ranging from 0 to 1
	# include the rownames and colnames into the matrix
	mat <- cbind(rownames(mat),mat)
	mat <- rbind(colnames(mat),mat)
	plot.window(xlim=c(1,ncol(mat)),ylim=c(1,nrow(mat)))
	# for line wise plotting
	switch(name,
			gamma=expr<-expression(Gamma),
			phi=expr<-expression(Phi),
			theta=expr<-expression(Theta),
			expr <- "Matrix")
	row.steps <- seq(1,nrow(mat))
	col.steps <- seq(1,ncol(mat))
	for(i in row.steps) { # rows
		for(j in col.steps) { # cols
			if(i==1 && j==1) {
				text(j,nrow(mat)-i+1,expr)
			} else {
				text(j,nrow(mat)-i+1,mat[i,j])	
			}
		}
	}
}

#construct time shifts in the gamma matrix
get.partitions <- function(TT,npairs) {
	partitions <- NULL
	for(i in 1:npairs) {
		shiftmax <- TT
		shifts <- NULL
		for(j in 1:TT) {
			if(shiftmax==0) {
				sh <- 0
			} else {
				sh <- sample(0:shiftmax,1,prob=c(0.000001,dexp(seq(1,shiftmax),0.6)))				
			}
			shifts <- c(shifts,sh)
			shiftmax <- shiftmax - sh
		}
		# permute to get all columns with equal probability shifted
		shifts <- sample(shifts)
		if(sum(shifts)!=TT) {
			i <- i-1
			next
		}
		if(is.null(partitions)) {
			partitions <- rbind(partitions, shifts)
		} else {
			if(length(which(apply(partitions,1,paste,collapse="_")==paste(shifts,collapse="_")))>0) {
				i <- i-1
			} else {
				partitions <- rbind(partitions, shifts)
			}
		}		
	}
	rownames(partitions) <- 1:nrow(partitions)
	colnames(partitions) <- 1:TT
	partitions
}


# for a given vector partition of the same column number as gammax, extract
# the columns as specified in partition:
# 10 0 0 0 0 0 0 0 0 would mean to take 10 times the first column
# 9  1 0 0 0 0 0 0 0 take 9 times the first column, 1 time the second column,
# etc.
get.gammashifted <- function(gammax, partitions) {
	gammax <- uniquegamma(gammax)
	stopifnot(ncol(gammax)==length(partitions))
	gammashift <- NULL
	for(i in 1:ncol(gammax)) {
		if(partitions[i]!=0) {
			gammashift <- cbind(gammashift,matrix(rep(gammax[,i],partitions[i]),ncol=partitions[i]))
		}
	}
	dimnames(gammashift) <- dimnames(gammax)
	gammashift
}

# written by Holger Froehlich, published in package 'nem' available on bioconductor
transitive.reduction <- function (g) 
{
	if (!(class(g) %in% c("matrix", "graphNEL"))) 
		stop("Input must be an adjacency matrix or graphNEL object")
	if (class(g) == "graphNEL") {
		g = as(g, "matrix")
	}
	g = g - diag(diag(g))
	type = (g > 1) * 1 - (g < 0) * 1
	for (y in 1:nrow(g)) {
		for (x in 1:nrow(g)) {
			if (g[x, y] != 0) {
				for (j in 1:nrow(g)) {
					if ((g[y, j] != 0) && (g[x, j] != 0) & (sign(type[x, 
												j]) * sign(type[x, y]) * sign(type[y, j]) != 
								-1)) 
						g[x, j] = 0
				}
			}
		}
	}
	g
}
# written by Holger Froehlich, published in package 'nem' available on bioconductor
transitive.closure <- function (g, mat = FALSE, loops = TRUE) 
{
	if (!(class(g) %in% c("graphNEL", "matrix"))) 
		stop("Input must be either graphNEL object or adjacency matrix")
	g <- as(g, "matrix")
	n <- ncol(g)
	matExpIterativ <- function(x, pow, y = x, z = x, i = 1) {
		while (i < pow) {
			z <- z %*% x
			y <- y + z
			i <- i + 1
		}
		return(y)
	}
	h <- matExpIterativ(g, n)
	h <- (h > 0) * 1
	dimnames(h) <- dimnames(g)
	if (!loops) 
		diag(h) <- rep(0, n)
	else diag(h) <- rep(1, n)
	if (!mat) 
		h <- as(h, "graphNEL")
	return(h)
}
# for timepoints that are replicated
# timepoints in columns
tp.median <- function(dat) {
	levs <- unique(colnames(dat))
	dat.m <- NULL
	for(l in levs) {
		dd <- dat[,which(colnames(dat)==l)]
		dd.m <- apply(dd,1,median)
		dat.m <- cbind(dat.m,dd.m)
	}
	colnames(dat.m) <- levs
	return(dat.m)
} 

# transform multiple edge types to only one edge type
detailed.to.simple.regulations <- function(phi) {
	ifelse(phi==0,0,1)
	#phi[phi==2] <- 1
	#phi
}
