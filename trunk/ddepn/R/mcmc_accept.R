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
	#laplace <- !is.null(lambda) && !is.null(B) && !is.null(Z)
	#scalefree <- !is.null(gam) && !is.null(it) && !is.null(K)
	scalefac <- 1# 3*length(bettermodels[[1]]$phi)
##	pr <- bestmodel$pr
##	pnew <- bettermodels[[1]]$posterior / scalefac # scale by dimension of metropolis jumping distribution (hope this is right)
##	pold <- bestmodel$posterior / scalefac
##	diffpost <- pnew - pold
##	diffproposal <- bettermodels[[1]]$pegmundo - bettermodels[[1]]$pegm
#browser()
##	# scale down the difference to get reasonable acceptance rates
##	diffpost <- sign(diffpost) * log10(abs(diffpost))
##	if(is.na(diffpost))# if equal posteriors, then don't take
##		diffpost <- -1
##	diffproposal <- sign(diffproposal) * log10(abs(diffproposal))
##	if(is.na(diffproposal))
##		diffproposal <- -1
##	diffs <- diffpost + diffproposal #pnew-pold + bettermodels[[1]]$pegmundo-bettermodels[[1]]$pegm
##	acpt <- min(exp(diffs),1)
	
	pnew <- bettermodels[[1]]$posterior
	prnew <- bettermodels[[1]]$pr
	pold <- bestmodel$posterior
	prold <- bestmodel$pr
	postratio <- pnew/pold
	prratio <- prnew/prold
	lratio <- bettermodels[[1]]$L / bestmodel$L
	#p2 <- prnew/prold
#browser()
	#diffproposal <- bettermodels[[1]]$pegmundo/bettermodels[[1]]$pegm
	postratio <- pnew/pold
	diffproposal <- bettermodels[[1]]$pegmundo / bettermodels[[1]]$pegm
	#postratio <- exp(pnew - pold)
	#diffproposal <- exp(bettermodels[[1]]$pegmundo - bettermodels[[1]]$pegm)
	acpt <- min(postratio*diffproposal,1)
	
	if(acpt==Inf || bettermodels[[1]]$posterior==-Inf || bettermodels[[1]]$posterior==Inf) {
		browser()
	}
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
				revert=type.txt<-"reverting edge")
		if(bestmodel$posterior < bestproposal$posterior) {
			#m <- paste("+++", type.txt, " MAP old: ", signif(pold,digits=7), " MAP new: ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
			m <- paste(type.txt, " MAP old/new: ", signif(pold,digits=7), " / ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
		} else {
			#m <- paste("~~~", type.txt, " MAP old: ", signif(pold,digits=7), " MAP new: ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
			m <- paste(type.txt, " MAP old/new: ", signif(pold,digits=7), " / ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
		}
	} else {
		bestproposal <- bestmodel
		#m <- paste("---MAP old: ", signif(pold,digits=7), " MAP new: ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
		m <- paste("MAP old/new: ", signif(pold,digits=7), " / ", signif(pnew,digits=7), " acpt: ", round(acpt,digits=5))
	}
	## acceptance newlambda
	if(priortype=="laplaceinhib" || priortype=="laplace") {
		pnew <- posterior(bestproposal$phi, bestproposal$L, newlambda, B, Z, gam, it, K, priortype) / scalefac	
	} else if(priortype=="scalefree") {
		pnew <- posterior(bestproposal$phi, bestproposal$L, NULL, B, Z, newlambda, it, K, priortype) / scalefac
	}
	
	postratio2 <- pnew/pold
	##diffpost2 <- pnew - pold
	####diffpost2 <- sign(diffpost2) * log10(abs(diffpost2))
	####if(is.na(diffpost2))# if equal posteriors, than take with p=0.5
	####	diffpost2 <- -1	
	##lacpt <- min(exp(diffpost2),1)
	lacpt <- min(postratio2, 1)
	m <- paste(m, " lacpt: ", round(lacpt,digits=5), "MAP ratio:",signif(postratio,digits=3), "Prior ratio:", signif(prratio,digits=3))

	takeit <- sample(c(0,1),1,prob=c((1-lacpt),lacpt))
	if(takeit==1) {
		if(priortype=="laplaceinhib" || priortype=="laplace") {
			bestproposal$lambda <- newlambda
		} else if(priortype=="scalefree") {
			bestproposal$gam <- newlambda
		}
		bestproposal$posterior <- pnew #newpost
	}
	if(priortype=="laplaceinhib" || priortype=="laplace") {
		m <- paste(m, " lambda: ", bestproposal$lambda)
	} else if(priortype=="scalefree") {
		m <- paste(m, " gam: ", bestproposal$gam)
	}
	
	print(m)
	return(list(bestproposal=bestproposal, acpt=acpt, lacpt=lacpt))
}


