# Given a gamma state transition matrix, expand to a number of replicates.
# 
# Author: benderc
###############################################################################

expand.gamma <- function(gammaposs, reps=4, order=TRUE) {
	gamma.ex <- NULL
	for(i in 1:reps) {
		gamma.ex <- cbind(gamma.ex, gammaposs)
	}
	if(order) {
		ord <- order(as.numeric(colnames(gamma.ex)))
		gamma.ex <- gamma.ex[,ord]		
	}
	gamma.ex
}

