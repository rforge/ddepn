# reformat colnames of example dataset to conform to ddepn input matrix requirements
# 
# Author: benderc
###############################################################################

format_ddepn <- function(dat) {
	colnames(dat) <- sub("-[0-9]+_","_",colnames(dat))
	dat
}

