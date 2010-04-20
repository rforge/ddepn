# Generation of artificial data.
# needs package KEGGgraph
# Author: benderc
###############################################################################


# given a network phi, simulate artificial timecourse data
simulatedata <- function(phi, mu.bg=0, sd.bg=0.1,
		mu.signal.a=1, sd.signal.a=0.5,
		mu.signal.i=-1, sd.signal.i=0.5,
		stimulus=sample(nrow(phi),2),TT=10,R.t=4,R.b=3,
		plot=FALSE){
	print(paste("Stimuli: ",paste(stimulus,rownames(phi)[stimulus],sep="->"),sep=" "))
	tps <- 1:TT
	reps <- R.t*R.b
	gammaposs <- uniquegammaposs(propagate.effect.simple(phi,1:TT,stimulus,reps))
	gammax <- gammaposs[,sort(sample(1:ncol(gammaposs), TT, replace=T))]
	colnames(gammax) <- tps
	gammax <- expand.gamma(gammax,reps)
	DD <- get.data(gammax,mu.bg,sd.bg,mu.signal.a,sd.signal.a,mu.signal.i,sd.signal.i,stimulus,TT,R.t,R.b)
	datx <- DD$datx
	downreg <- DD$downreg
	
	if(plot) {
		datx.m <- tp.median(datx)
		par(mfrow=c(2,2))
		plotdetailed(phi)
		plotmatrix(phi,"phi")
		plotmatrix(uniquegamma(gammax),"gamma")
		cols <- rainbow(nrow(datx.m))
		leg <- rownames(datx.m)
		ylim <- range(datx.m)
		ylim[2] <- ylim[2] +  2*abs(ylim[2])
		plot(colnames(datx.m),datx.m[1,],type='l',col=cols[1],ylim=ylim)
		sapply(1:nrow(datx.m), function(i,datx.m,cols) lines(colnames(datx.m),datx.m[i,],col=cols[i]), datx.m=datx.m,cols=cols)
		legend("topleft",legend=leg,fill=cols,lty=1,cex=0.7,ncol=2)
	}
	return(list(datx=datx,phix=phi,gammax=gammax,stimx=stimulus,tps=tps,R.t=R.t,R.b=R.b,downreg=downreg, gammaposs=gammaposs))
}


