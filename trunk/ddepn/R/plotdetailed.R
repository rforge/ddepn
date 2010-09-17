# Plot graphs containing inhibitions
# input is an adjacency matrix phi with
#  phi[i,j]==1 : activation
#  phi[i,j]==2 : inhibition
#  phi[i,j]==0 : no edge from i to j
# 
# Author: benderc
###############################################################################

get.labels <- function(phi) {
	labels <- NULL
	l.names <- NULL
	for(i in 1:nrow(phi)) {
		for(j in 1:ncol(phi)) {
			#if(phi[i,j]==0) {
			#	next
			#}
			labels <- c(labels,phi[i,j])
			l.names <- c(l.names, paste(rownames(phi)[i],colnames(phi)[j],sep="~"))
		}
	}
	names(labels) <- l.names
	labels
}
get.arrowhead <- function(phi) {
	arrowhead <- NULL
	l.names <- NULL	
	for(i in 1:nrow(phi)) {
		for(j in 1:ncol(phi)) {
			#if(phi[i,j]==0) {
			#	next
			#}
			if(phi[i,j]==3) {
				# can be both activation or inhibition
				arrowhead <- c(arrowhead,"odot")				
			} else {
				if(phi[i,j]==2){
					arrowhead <- c(arrowhead,"tee")
				} else {
					arrowhead <- c(arrowhead,"open")	
				}
			}
			l.names <- c(l.names, paste(rownames(phi)[i],colnames(phi)[j],sep="~"))
		}
	}
	names(arrowhead) <- l.names
	arrowhead
}
get.arrowtail <- function(phi) {
	arrowtail <- NULL
	l.names <- NULL
	for(i in 1:nrow(phi)) {
		for(j in 1:ncol(phi)) {
			arrowtail <- c(arrowtail,"none")
			l.names <- c(l.names, paste(rownames(phi)[i],colnames(phi)[j],sep="~"))
		}
	}
	names(arrowtail) <- l.names
	arrowtail
}
plotdetailed <- function(phi, weights=NULL,main="",stimuli=NULL, layoutType="dot", fontsize=20) {
	phix <- phi
	if(all(phix==0)) {
		phix[1,1] <- 1
		emptygraph <- TRUE
	} else {
		emptygraph <- FALSE
	}
	phix[which(phix==2)] <- 1
	g1 <- as(phix,"graphNEL")
	if(!is.null(weights)) {
		labels <- get.labels(weights)
	} 
	arrowhead <- get.arrowhead(phi)
	graph.par(list(graph = list(main = main)), edges = list(lwd = 1)) #, sub = "... and a subtitle",cex.main = 1.8, cex.sub = 1.4, col.sub = "gray")))
	if(emptygraph)
		graph.par(list(edges = list(lwd = 0)))
	## something with the edgeRenderInfo function doesn't work, so pass edgeAttrs directly
	## to layoutGraph
	## In stead, passing nodeAttrs directly to layoutGraph doesn't seem to work, so
	## use the nodeRenderInfo function here
	#edgeRenderInfo(g1) <- list(label = labels, arrowhead=arrowhead) # , arrowtail=arrowtail)
	if(is.null(stimuli)) {
		nodeRenderInfo(g1) <- list(shape = "box", fill="lightgray", fontsize=fontsize)
	} else {
		fills <- rep("lightgray",length(nodes(g1)))
		names(fills) <- nodes(g1)
		fills[unique(names(unlist(stimuli)))] <- "red"
		nodeRenderInfo(g1) <- list(shape = "box", fill=fills, fontsize=fontsize)
	}
	if(!is.null(weights)) {
		if(emptygraph) {
			g1 <- layoutGraph(g1, layoutType=layoutType) #, recipEdges = "distinct", layoutType=layoutType)
		} else {
			g1 <- layoutGraph(g1, recipEdges = "distinct", edgeAttrs = list(label = labels, 
							arrowhead = arrowhead), layoutType=layoutType)
		}
	} else {
		g1 <- layoutGraph(g1, recipEdges = "distinct", edgeAttrs = list(arrowhead = arrowhead),
						layoutType=layoutType)
	}
	renderGraph(g1)
	if(emptygraph)
		graph.par(list(edges = list(lwd = 1)))
}
