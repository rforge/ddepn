# TODO: Add comment
# 
# Author: benderc
###############################################################################
#
# wrapper for substitute
#
esub <- function(expr, sublist) {
	do.call("substitute", list(expr, sublist))
}

bestgam <- function(y, tp, xn, selection.criterion="aic") {
	## get the fits and store them in a list
	numPoints <- length(unique(tp))-2
	gamres=vector("list", length=numPoints+1)
	gamres[[1]] = gam(y~1) # constant fit
	gamres[[2]] = gam(y~tp) # straight line fit
	for(i in 2:length(gamres)-1) {	
		# create an expression for the gam function
		# first use quote to assemble the expression
		expre <- quote(y~s(tp,i))
		# use the substitute wrapper to change the i into its value
		# as.numeric is needed, otherwise eval returns an integer with the "L" notation
		# which gives an error in the gam function
		expre <- esub(expre, list(i=eval(as.numeric(i))))
		gamres[[i+1]] = gam(eval(expre))
	}
	## color vector for the different fits
	ccc <- gray(0:length(gamres) / length(gamres))
	## test which fit is the best
	aovtab <- NULL
	for(model in 1:length(gamres)) {
		aovtab <- rbind(aovtab,(as.matrix(anova(gamres[[1]],gamres[[model]],test="Chisq"))[2,]))
	}
	# use maximal aic or bic as model selection criterion
	pval <- aovtab[-1,5]
	df <- aovtab[-1,3]
	N <- length(y)
	if(selection.criterion=="aic") {
		aic <- 2*df - 2*log(pval) # maximise aic
		i <- min(which(aic==max(aic)))
	} else if(selection.criterion=="aic") {
		bic <- -2*log(pval) + df * log(N)
		i <- min(which(bic==max(bic))) # maximise bic
	} else {
		i <- max(which(pval==min(pval))) # minimise p-value
	}
	p <- aovtab[(i+1),5]
	main <- ""
	if(p<=0.05) {
		main <- c(main, paste("Best fit: Spline df=", i, " p=", signif(p,3), sep=""))
	} else {
		i <- 1
		main <- c(main,paste("No significant change, p=",signif(p,3)))
	}
	pred <- predict(gamres[[i]], newdata=data.frame(tp=xn))
	return(pred)
}

