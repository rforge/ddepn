# TODO: plot wrapper to give a summary plot after inference
# 
# Author: benderc
###############################################################################

plotrepresult <- function(lst,outfile=NULL,fontsize=15) {
	phi.final <- lst$phi
	phi.orig <- lst$phi.orig
	weights <- lst$weights
	stats <- lst$stats
	stimuli <- lst$stimuli
	dat <- lst$dat
	phi.final.tr <- transitive.reduction(phi.final)
	if(!is.null(outfile)) {
		pdf(outfile)
	}
	par(mfrow=c(1,3))
	if(all(phi.orig==0)) {
		plot.new()
		text(x=0.5,y=0.5,"No reference network given.")
	} else {
		plotdetailed(phi.orig,main="Original Graph",stimuli=stimuli,fontsize=fontsize)		
	}
	plotdetailed(phi.final,weights,main="Inferred Graph",stimuli=stimuli,fontsize=fontsize)
	#plotdetailed(phi.final.tr,main="Inferred, transitive.reduction",stimuli=stimuli,fontsize=fontsize)	
	#boxplot(as.data.frame(t(dat)),las=2,main="Dataset's data distribution")
	if(all(phi.orig==0)) {
		plot.new()
		text(x=0.5,y=0.5,"No SN/SP plot possible.")
	} else {
		if(!is.null(lst$burnin)) {
			#sts <- as.data.frame(stats[-(1:(lst$burnin*4)),c("sn","sp")])
			s1 <- min(which(stats[,"sp"]!=1 & stats[,"sp"]!=0))
			s2 <- min(which(stats[,"sn"]!=0 & stats[,"sn"]!=1))
			st <- max(s1,s2)
			st <- min(st,lst$burnin)
			sts <- as.data.frame(stats[-(1:st),c("sn","sp")])
		} else {
			sts <- as.data.frame(stats[,c("sn","sp")])
		}
		boxplot(sts,ylim=c(0,1),main="SN/SP")
		legend("bottomleft",legend=c("sn=tp/tp+fn", "sp=tn/tn+fp"))
	}
	if(!is.null(outfile)) {
		dev.off()
	}
}
