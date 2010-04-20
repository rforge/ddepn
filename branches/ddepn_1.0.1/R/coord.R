# given a matrix index k, find the corresponding coordinate in the matrix mat
# 
# Author: benderc
###############################################################################


coord <- function(k,mat) {
	y <- ceiling(k/nrow(mat))
	x <- k%%nrow(mat)
	if(x==0) {
		x <- nrow(mat)
	}
	return(c(x,y))
}

