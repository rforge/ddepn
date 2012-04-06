# Test Suite for the ddepn package
#
# Author: benderc
###############################################################################

## clear the console
cls <- function() {	
	system("clear")
} 

wrapddepn <- function(dataset, phit, stimuli, debug=FALSE, plotresults=TRUE, ...) {
	burnin <- 5
	## if debug==TRUE, options(error=recover) can be used, but not if the ddepn call
	## is wrapped in the try statement
	if(debug) {
		msg <- ddepn(dataset$datx, phiorig=phit, maxiterations=20, p=15, q=0.3, m=0.8, burnin=burnin, plotresults=plotresults, ...)
	} else {
		msg <- try(ddepn(dataset$datx, phiorig=phit, maxiterations=20, p=15, q=0.3, m=0.8, , burnin=burnin, plotresults=plotresults, ...))
		if(class(msg)=="try-error") {
			msg <- "failed"
		} else {
			## extra treatment of mcmc testing, since return object 
			## is a list containing the samplings and traces
			if("samplings" %in% names(msg)) {
				print("inhibMCMC done.")
				if(plotresults)
					plot_mcmctraces(msg,thin=10)
				#browser()
			} else { # GA was used
				print("GA done.")
				#if(plotresults)
				#	plotdetailed(msg$phi,stimuli=stimuli,fontsize=20)
			}
			msg <- "passed"
		}
	}
	invisible(msg)
}

run_checks <- function(outfile=NULL, plotresults=TRUE) {
	## sample a network
	set.seed(1234)
	n <- 6
	signet <- signalnetwork(n=n, nstim=2, cstim=0, prop.inh=0.2)
	phit <- signet$phi
	stimuli <- signet$stimuli
	
	## sample data
	dataset <- makedata(phit, stimuli, mu.bg=1200, sd.bg=400, mu.signal.a=2000, sd.signal.a=1000)
	
	## use original network as prior matrix
	## reset all entries for inhibiting edges 
	## to -1
	B <- phit
	B[B==2] <- -1
	
	## plots will be stored in a file, ask for a name if not given
	## if empty input, then use default name
	if(plotresults) {
		if(is.null(outfile)) {
			defaultstr <- "testddepn.pdf"
			outfile <- readline(prompt = paste("Store output in file [",defaultstr,"]?",sep=""))
			if(outfile=="") {
				outfile <- defaultstr
			}
		}
	}
	tests <- c("1.1. GA, BIC optimisation",
			"1.2. GA, uniform prior",
			"1.3.1 GA, laplaceinhib, two types, lambda fix",
			"1.3.2 GA, laplaceinhib, two types, lambda integrated",
			"1.4.1 GA, laplace, one type, lambda fix",
			"1.4.2 GA, laplace, one type, lambda integrated",
			"1.5. GA, scalefree",
			"2.1. inhibMCMC, uniform prior",
			"2.2.1 inhibMCMC, laplaceinhib, two types, lambda fix",
			"2.2.2 inhibMCMC, laplaceinhib, two types, lambda integrated",
			"2.3.1 inhibMCMC, laplace, one type, lambda fix",
			"2.3.2 inhibMCMC, laplace, one type, lambda integrated",
			"2.4. inhibMCMC, scalefree",
			"2.5. inhibMCMC, laplaceinhib, sample lambda")
	msgs <- NULL
	
	## open the pdf connection
	if(plotresults)
		pdf(outfile)
	############
	## Genetic algorithm, 
	
	## using BIC score as optimisation criterion
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[1]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=TRUE, plotresults=plotresults))
	
	## using a uniform prior
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[2]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=FALSE, priortype="uniform", plotresults=plotresults))
	
	## using laplaceinhib, lambda fix
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[3]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=FALSE, priortype="laplaceinhib", lambda=0.01, B=B, plotresults=plotresults))
	
	## using laplaceinhib, lambda integrated
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[4]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=FALSE, priortype="laplaceinhib", lambda=NA, B=B, plotresults=plotresults))
	
	## using laplace, lambda fix
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[5]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=FALSE, priortype="laplace", lambda=0.01, B=B, plotresults=plotresults))
	
	## using laplace, lambda integrated
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[6]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=FALSE, priortype="laplace", lambda=NA, B=B, plotresults=plotresults))
	
	## using laplace, lambda integrated
	cls()
	print(paste("1. Running GA tests."))
	print(paste("Running", tests[7]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="netga", usebics=FALSE, priortype="scalefree", gam=2.2, it=500, K=0.8, plotresults=plotresults))
	
	
	
	############
	## MCMC sampling using a uniform prior
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[8]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="uniform", plotresults=plotresults))
	
	## laplaceinhib
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[9]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="laplaceinhib", lambda=0.01, B=B, plotresults=plotresults))
	
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[10]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="laplaceinhib", lambda=NA, B=B, plotresults=plotresults))
	
	## laplace
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[11]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="laplace", lambda=0.01, B=B, plotresults=plotresults))
	
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[12]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="laplace", lambda=NA, B=B, plotresults=plotresults))
		
	## MCMC sampling using a scale free prior
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[13]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="scalefree", gam=2.2, it=500, K=0.8, plotresults=plotresults))
	
	## laplace inhib, sample the lambda parameter
	cls()
	print(paste("2. Running MCMC tests."))
	print(paste("Running", tests[14]))
	msgs <- c(msgs, wrapddepn(dataset, phit, stimuli, inference="mcmc", usebics=FALSE, priortype="laplaceinhib", lambda=0.01, samplelambda=0.1, B=B, plotresults=plotresults))
	
	## close pdf
	if(plotresults)
		dev.off()
		
	report <- cbind(tests, msgs)
	print(report)
	invisible(report)
}
