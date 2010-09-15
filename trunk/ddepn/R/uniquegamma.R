# Finding unique state-vector matrices
# 
# Author: benderc
###############################################################################


## if many replicates for each timepoint, the gamma matrix will have multiple equal columns
## -> reduce to the unique gamma matrix
uniquegamma <- function(gammax) {
	#ind.uniquegamma <- 1 + (as.numeric(unique(colnames(gammax)))-1) * length(colnames(gammax))/length(unique(colnames(gammax)))
	reps <- length(colnames(gammax))/length(unique(colnames(gammax)))
	ind.uniquegamma <- seq(1,ncol(gammax),by=reps)
	gammax[,ind.uniquegamma]
}

## for multiple experiments, the gammaposs-matrix contains all possible state vectors
## reduce the matrix with duplicate columns such that the order of the states is preserved
uniquegammaposs <- function(gammaposs.tmp) {
	gammaposs <- gammaposs.tmp[,1,drop=F]
	cols <- colnames(gammaposs.tmp)[1]
	if(ncol(gammaposs.tmp)>1) {
		for(i in 2:ncol(gammaposs.tmp)) {
			if(all(gammaposs.tmp[,i]==gammaposs.tmp[,i-1]))
				next
			gammaposs <- cbind(gammaposs,gammaposs.tmp[,i])
			cols <- c(cols, colnames(gammaposs.tmp)[i])
		}
	}
	colnames(gammaposs) <- cols
	return(gammaposs)
}

