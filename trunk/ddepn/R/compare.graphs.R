# Count number of correctly/falsly inferred edges
# 
# Author: benderc
###############################################################################


compare.graphs.tc <- function(phiorig,phi,ignore.type=TRUE) {
	if(ignore.type) {
		originalnc <- detailed.to.simple.regulations(phiorig)
		matnc <- detailed.to.simple.regulations(phi)
	} else {
		originalnc <- phiorig
		matnc <- phi
	}
	diag(matnc) <- 0
	diag(originalnc) <- 0

	## inhibitionen/aktivierungen richtig z�hlen
	tp <- length(which(originalnc==matnc & (originalnc==1 | originalnc==2)))
	tn <- length(which(originalnc==matnc & originalnc==0))
	fn <- length(which(originalnc!=matnc & (originalnc==1 | originalnc==2)))
	fp <- length(which(originalnc!=matnc & originalnc==0))
		
	tpr = tp/(tp+fn) # sn
	tnr = tn/(tn+fp) # sp
	prec <- round(tp/(tp+fp),digits=4)
	f1 <- round(2*prec*tpr/(prec+tpr),digits=4)
	return(data.frame(tp=tp,tn=tn,fp=fp,fn=fn,sn=tpr,sp=tnr,prec=prec,f1=f1))	
}
