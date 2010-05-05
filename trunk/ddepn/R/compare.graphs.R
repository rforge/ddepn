# Count number of correctly/falsly inferred edges
# 
# Author: benderc
###############################################################################


compare.graphs.tc <- function(O,M) {
	#originalnc <- detailed.to.simple.regulations(O)
	#originalc <- transitive.closure(detailed.to.simple.regulations(O),mat=T)
	#matnc <- detailed.to.simple.regulations(M)
	#matc <- transitive.closure(detailed.to.simple.regulations(M),mat=T)
	originalnc <- O
	matnc <- M
	#diag(matc) <- 0
	#diag(originalc) <- 0
	diag(matnc) <- 0
	diag(originalnc) <- 0

	## inhibitionen/aktivierungen richtig zählen
	tp <- length(which(originalnc==matnc & (originalnc==1 | originalnc==2)))
	tn <- length(which(originalnc==matnc & originalnc==0))
	fn <- length(which(originalnc!=matnc & (originalnc==1 | originalnc==2)))
	fp <- length(which(originalnc!=matnc & originalnc==0))
	
	## sieht gut aus, nicht sicher, ob ganz richtig, addiert sich jedenfalls zu
	## tp + tn + fp + fn = summe(kanten)
	## gründe für dieses zählen:
	## tp: kanten, die implizit in O vorhanden, werden oft direkt in I gelernt, daher tp
	## fp: kanten, die direkt in I gelernt werden, nur dann falsch, wenn auch nicht implizit in O
	## tn: wenn methode sagt, in den daten gibt es die kanten nicht (Iij=0) UND die kante ist
	##     in dem originalgraphen auch wirklich nicht da, dann ist das als TN zu zählen
	## fn: wenn ich kante in I nicht lerne, sie aber im originalgraphen drin ist, dann ist sie fn
	#tp = sum(matnc[originalc==1]==1,na.rm=TRUE)
	#tn = sum(matnc[originalnc==0]==0,na.rm=TRUE)
	#fp = sum(matnc[originalc==0]==1,na.rm=TRUE)
	#fn = sum(matnc[originalnc==1]==0,na.rm=TRUE) 

	## graphen soi lassen wie sie sind und nur echte, direkte verbindungen zählen
	## ist nicht besonders gut, aber immerhin besser als AUC 0.5
	#tp = sum(matnc[originalnc==1]==1,na.rm=TRUE)
	#tn = sum(matnc[originalnc==0]==0,na.rm=TRUE) 
	#fp = sum(matnc[originalnc==0]==1,na.rm=TRUE) 
	#fn = sum(matnc[originalnc==1]==0,na.rm=TRUE)  
	
	## wäre auch richtig:
	#tp = sum(matnc[originalc==1]==1,na.rm=TRUE)
	#tn = sum(matnc[originalc==0]==0,na.rm=TRUE) 
	#fp = sum(matnc[originalc==0]==1,na.rm=TRUE) 
	#fn = sum(matnc[originalc==1]==0,na.rm=TRUE)  

	# von holger vorgeschlagen,  addiert sich nicht auf summe(kanten)
	#tp = sum(matnc[originalnc==1]==1,na.rm=TRUE)
	#tn = sum(matnc[originalc==0]==0,na.rm=TRUE)
	#fp = sum(matnc[originalc==0]==1,na.rm=TRUE)
	#fn = sum(matnc[originalnc==1]==0,na.rm=TRUE)	
	
	tpr = tp/(tp+fn) # sn
	tnr = tn/(tn+fp) # sp
	prec <- round(tp/(tp+fp),digits=4)
	f1 <- round(2*prec*tpr/(prec+tpr),digits=4)
	return(data.frame(tp=tp,tn=tn,fp=fp,fn=fn,sn=tpr,sp=tnr,prec=prec,f1=f1))	
}
## ist kacke, gibt sehr hohe sn und niedrige sp 
#	tp = sum(matc[originalc==1]==1,na.rm=TRUE)
#	tn = sum(matc[originalc==0]==0,na.rm=TRUE) 
#	fp = sum(matc[originalc==0]==1,na.rm=TRUE) 
#	fn = sum(matc[originalc==1]==0,na.rm=TRUE)  


#
#
### zähle echte/falsche aktiv/passiv stati
#compare.graphs.ep <- function(O,M,tps,stimuli) {
#	oprop <- propagate.effect.set(O,tps,stimuli,reps=1)
#	mprop <- propagate.effect.set(M,tps,stimuli,reps=1)
#	tp = sum(mprop[oprop==1]==1,na.rm=TRUE)
#	tn = sum(mprop[oprop==0]==0,na.rm=TRUE) 
#	fp = sum(mprop[oprop==0]==1,na.rm=TRUE) 
#	fn = sum(mprop[oprop==1]==0,na.rm=TRUE)  
#	tpr = tp/(tp+fn) # sn
#	tnr = tn/(tn+fp) # sp
#	prec <- round(tp/(tp+fp),digits=4)
#	f1 <- round(2*prec*tpr/(prec+tpr),digits=4)
#	return(data.frame(tp=tp,tn=tn,fp=fp,fn=fn,sn=tpr,sp=tnr,prec=prec,f1=f1))	
#}
#	tp = sum(matnc[originalnc==1]==1,na.rm=TRUE)
#	tn = sum(matnc[originalc==0]==0,na.rm=TRUE) ## implizite kante in referenz soll nicht als negative kanten angesehen werden
#	fp = sum(matnc[originalc==0]==1,na.rm=TRUE) ## implizite kanten in referenz soll icht als negative kante angesehen werden
#	fn = sum(matnc[originalnc==1]==0,na.rm=TRUE) ##
##
#	tp = sum(matnc[originalc==1]==1,na.rm=TRUE)
#	tn = sum(matnc[originalnc==0]==0,na.rm=TRUE) ## implizite kante in referenz soll nicht als negative kanten angesehen werden
#	fp = sum(matnc[originalnc==0]==1,na.rm=TRUE) ## implizite kanten in referenz soll icht als negative kante angesehen werden
#	fn = sum(matnc[originalc==1]==0,na.rm=TRUE) ## 
##

#compare.graphs2 <- function(original,mat) {
#	origclosed <- transitive.closure(original,mat=T)
#	diag(mat) <- 0
#	diag(origclosed) <- 0
#	tp = sum(mat[original==1]==1,na.rm=TRUE)
#	tn = sum(mat[origclosed==0]==0,na.rm=TRUE)
#	fp = sum(mat[origclosed==0]==1,na.rm=TRUE)
#	fn = sum(mat[original==1]==0,na.rm=TRUE)	
#	tpr = tp/(tp+fn) # sn
#	tnr = tn/(tn+fp) # sp
#	prec <- round(tp/(tp+fp),digits=4)
#	f1 <- round(2*prec*tpr/(prec+tpr),digits=4)
#	return(data.frame(tp=tp,tn=tn,fp=fp,fn=fn,sn=tpr,sp=tnr,prec=prec,f1=f1))
#}
