# take each possible move as proposition and sample from these
# 
# Author: benderc
###############################################################################

random.gammaprime <- function(gammax) {
	gammax.u <- uniquegamma(gammax)
	cols <- sample(1:ncol(gammax.u),ncol(gammax.u),replace=T)
	gammax.u.rnd <- gammax.u[,cols]
	ret <- expand.gamma(gammax.u.rnd[,order(as.numeric(colnames(gammax.u.rnd)))],reps=ncol(gammax)/ncol(gammax.u))
	colnames(ret) <- colnames(gammax)
	return(ret)
}

