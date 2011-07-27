# TODO: Add comment
# 
# Author: benderc
###############################################################################


makedata <- function(phi, stimuli, R.t=3, R.b=3, TT=10,
		mu.bg=0, sd.bg=0.1,
		mu.signal.a=2, sd.signal.a=0.5,
		mu.signal.i=-2, sd.signal.i=0.5, allow.stim.off=FALSE) {
	dat <- gam <- NULL
	for(i in 1:length(stimuli)) {
		stimulus <- stimuli[[i]]
		datx <- simulatedata(phi, mu.bg=mu.bg, sd.bg=sd.bg,
				mu.signal.a=mu.signal.a, sd.signal.a=sd.signal.a,
				mu.signal.i=mu.signal.i, sd.signal.i=sd.signal.i,
				stimulus=stimulus,TT=TT,R.t=R.t,R.b=R.b,plot=F,stimuli=stimuli,allow.stim.off=allow.stim.off)
		gamt <- datx$gammax
		datt <- datx$datx
		for(j in unique(unlist(stimuli))) {
			datt[j,] <- rep(0,ncol(datt))	
		}
		sts <- paste(names(stimuli[[i]]),collapse="&")
		colnames(datt) <- colnames(gamt) <-  paste(sts,colnames(datt),sep="_")
		dat <- cbind(dat, datt)
		gam <- cbind(gam,gamt)	
	}
	dataset <- list(datx=dat, gammax=gam, stimuli=stimuli, phi=phi, R.t=R.t, R.b=R.b, TT=TT)
	dataset
}

