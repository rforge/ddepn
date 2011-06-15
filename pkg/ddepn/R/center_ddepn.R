# perform centering of a data matrix:
# subtract median of each time course (i.e. for each technical replicate) and 
# add the median of all biological and technical replicates of one experiment
#
# Author: benderc
###############################################################################


center_ddepn <- function(dat) {
	biolrepliid <- sub("_[0-9]*$","",colnames(dat))
	biolrepli <- unique(biolrepliid)
	experimentsid <- sub("-[0-9]*$","",biolrepliid)
	experiments <- unique(experimentsid)
	
	datc <- NULL
	for(exper in experiments) {		
		ind <- which(experimentsid==exper)
		xx <- dat[,ind]
		xxc <- NULL
		for(i in 1:nrow(xx)) {
			y <- xx[i,]
			tp <- sub("^.*_","",names(y))
			time <- unique(tp)
			brep <- as.numeric(sub("(.*-)([0-9]*)(_.*$)","\\2",names(y)))
			brepf <- unique(brep)
			mat <- matrix(y, length(tp)/length(time), length(time))
			mat <- mat - apply(mat, 1, median)
			backshift <- max(abs(1-min(mat)),median(y))
			mat <- mat + backshift
			if (any(as.vector(mat) <= 0)) {
				warning("Negative values occured during centering.")
			}
			yprime <- as.vector(mat)
			names(yprime) <- names(y)
			xxc <- rbind(xxc,yprime)
		}
		datc <- cbind(datc, xxc)
	}
	rownames(datc) <- rownames(dat)
	datc
}


