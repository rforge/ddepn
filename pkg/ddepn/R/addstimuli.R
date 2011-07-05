# TODO: Add comment
# 
# Author: benderc
###############################################################################

## add the stimuli to the data matrix
addstimuli <- function(dat) {
	cols <- colnames(dat)
	## first get the stimuli from the colnames
	tmp <- sapply(cols, function(x) strsplit(x, "_")[[1]])
	exps <- unique(tmp[1, ])
	allstim <- as.character(unique(unlist(sapply(exps, function(x) strsplit(x, "&")[[1]]))))
	stimuli <- list()
	for (i in 1:length(exps)) {
		expsi <- exps[i]
		stims <- as.character(sapply(expsi, function(x) strsplit(x, 
									"&")[[1]]))
		stimsid <- NULL
		for (j in 1:length(stims)) {
			el <- stims[j]
			x <- match(el, rownames(dat))
			if (any(is.na(x))) 
				x <- match(el, allstim)
			names(x) <- el
			stimsid <- c(stimsid, x)
		}
		names(stimsid) <- stims
		stimuli[[i]] <- stimsid
	}
	stimm <- match(unique(names(unlist(stimuli))), rownames(dat))
	## check if the stimuli nodes are present in the data matrix
	## if not, add them as dummy nodes.
	if (any(is.na(stimm))) {
		xx <- unlist(stimuli)
		xxmat <- unique(cbind(xx, names(xx)))
		toattach <- matrix(0, nrow = nrow(xxmat), ncol = ncol(dat), 
				dimnames = list(xxmat[, 2], colnames(dat)))
		prune <- stimm[!is.na(stimm)]
		if (length(prune) > 0) 
			dat <- dat[-prune, ]
		dat <- rbind(toattach, dat)
	}
	list(dat=dat,stimuli=stimuli)
}
