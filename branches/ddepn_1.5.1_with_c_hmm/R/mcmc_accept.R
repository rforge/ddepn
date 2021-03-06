# TODO: Add comment
# 
# Author: benderc
###############################################################################


mcmc_accept <- function(bestmodel, bettermodels, newlambda) {
	lambda <- bestmodel$lambda
	B <- bestmodel$B
	Z <- bestmodel$Z
	gam <- bestmodel$gam
	it <- bestmodel$it
	K <- bestmodel$K
	priortype <- bestmodel$priortype
	scalefac <- bestmodel$scalefac #0.005
	## get the scores
	pnew <- bettermodels[[1]]$posterior
	prnew <- bettermodels[[1]]$pr
	pold <- bestmodel$posterior
	prold <- bestmodel$pr
	
	## now the ratios (are differences since all scores are on log-scale
	prratio <- prnew - prold
	lratio <- bettermodels[[1]]$L - bestmodel$L  
	postratio <- pnew - pold
	postratio <- prratio + lratio
	
	## get the proposal distribution ratio
	diffproposal <- bettermodels[[1]]$pegmundo - bettermodels[[1]]$pegm
	
	## get acceptance rate, scaled by scalefactor that is inferred during the burnin
	acpt <- min(exp(scalefac * (postratio+diffproposal)),1)

	## did something go wrong??
	if(acpt==Inf || bettermodels[[1]]$posterior==-Inf || bettermodels[[1]]$posterior==Inf) {
		browser()
	}
	
	## accept or not
	takeit <- sample(c(0,1),1,prob=c((1-acpt),acpt))
	if(takeit==1) {
		bestproposal <- bettermodels[[1]]
		type <- bestproposal$lastmove
		switch(type,
				add=type.txt<-paste("adding edge"),
				addactivation=type.txt<-"adding activation",
				addinhibition=type.txt<-"adding inhibition",
				switchtype=type.txt<-"switching type",
				delete=type.txt<-"deleting edge",
				revert=type.txt<-"reverting edge",
				revswitch=type.txt<-"revswitch")
		## getting better
		if(bestmodel$posterior < bestproposal$posterior) {
			m <- paste("+++", type.txt, " MAP old/new: ", signif(pold,digits=7), " / ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
		} else { ## getting worse
			m <- paste(type.txt, " MAP old/new: ", signif(pold,digits=7), " / ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
		}
	} else {
		## if not taken, keep everything as it is.
		bestproposal <- bestmodel
		m <- paste("MAP old/new: ", signif(pold,digits=7), " / ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
	}
	## acceptance newlambda
	if(priortype %in% c("laplaceinhib","laplace","uniform")) {
		pnew <- posterior(bestproposal$phi, bestproposal$L, newlambda, B, Z, gam, it, K, priortype) #/ scalefac	
	} else if(priortype=="scalefree") {
		pnew <- posterior(bestproposal$phi, bestproposal$L, NULL, B, Z, newlambda, it, K, priortype) #/ scalefac
	}
	
	## posterior ratio for newlambda
	postratio2 <- pnew - pold
	## acceptance ratio for newlambda
	lacpt <- min(exp(scalefac*postratio2), 1)
	
	## update the output message
	m <- paste(m, " lacpt: ", round(lacpt,digits=5), "MAP ratio:",signif(postratio,digits=3), "Prior ratio:", signif(prratio,digits=3))

	## take newlambda or not
	takeit <- sample(c(0,1),1,prob=c((1-lacpt),lacpt))
	if(takeit==1) {
		if(priortype %in% c("laplaceinhib","laplace","uniform")) {
			bestproposal$lambda <- newlambda
		} else if(priortype=="scalefree") {
			bestproposal$gam <- newlambda
		}
		bestproposal$posterior <- pnew
	}
	if(priortype %in% c("laplaceinhib","laplace","uniform")) {
		m <- paste(m, " lambda: ", bestproposal$lambda)
	} else if(priortype=="scalefree") {
		m <- paste(m, " gam: ", bestproposal$gam)
	}
	m <- paste(m, " scale: ",scalefac,sep="")
	print(m)
	return(list(bestproposal=bestproposal, acpt=acpt, lacpt=lacpt))
}


