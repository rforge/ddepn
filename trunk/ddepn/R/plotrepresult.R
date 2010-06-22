# TODO: Add comment
# 
# Author: benderc
###############################################################################

plotrepresult <- function(lst,pdf=NULL,fontsize=15) {
	phi.final <- lst$phi
	phi.orig <- lst$phi.orig
	weights <- lst$weights
	result <- lst$result
	stimuli <- lst$stimuli
	dat <- lst$dat
	phi.final.tr <- transitive.reduction(phi.final)
	if(!is.null(pdf)) {
		pdf(pdf)
	}
	par(mfrow=c(2,2))
	if(all(phi.orig==0)) {
		plot.new()
		text(x=0.5,y=0.5,"No reference network given.")
	} else {
		plotdetailed(phi.orig,main="Original Graph",stimuli=stimuli,fontsize=fontsize)		
	}
	plotdetailed(phi.final,weights,main="Inferred Graph, weight>=0.8",stimuli=stimuli,fontsize=fontsize)
	plotdetailed(phi.final.tr,main="Inferred, transitive.reduction",stimuli=stimuli,fontsize=fontsize)	
	#boxplot(as.data.frame(t(dat)),las=2,main="Dataset's data distribution")
	if(all(phi.orig==0)) {
		plot.new()
		text(x=0.5,y=0.5,"No SN/SP plot possible.")
	} else {
		boxplot(as.data.frame(result[,c("sn","sp")]),ylim=c(0,1),main="SN/SP")
		legend("bottomleft",legend=c("sn=tp/tp+fn", "sp=tn/tn+fp"))
	}
	
#	if(!is.null(result)) {
#		if(is.null(pdf))
#			x11()
#		par(mfrow=c(1,3))
#		if(all(phi.orig==0)) {
#			plot.new()
#			text(x=0.5,y=0.5,"No SN/SP plot possible.")
#			plot.new()
#			text(x=0.5,y=0.5,"No TP/TN/FP/FN plot possible.")
#		} else {
#			boxplot(as.data.frame(result[,c("sn","sp")]),ylim=c(0,1))
#			legend("bottomleft",legend=c("sn=tp/tp+fn", "sp=tn/tn+fp"))
#			boxplot(as.data.frame(result[,c("tp","fp","tn","fn")]),ylim=c(0,max(result[,1:4])))
#		}
#		boxplot(as.data.frame(result[,c("elapsed")]),ylim=c(0,max(result[,c("elapsed")])),ylab="sec",xlab="elapsed time per run")	
#	}
	if(!is.null(pdf)) {
		dev.off()
	}
}


