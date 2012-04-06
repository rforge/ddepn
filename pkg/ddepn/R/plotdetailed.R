# source("~/projects/rfunctions/graph/plotdetailed_igraph.R")
# 
# Author: benderc
###############################################################################

layout.ellipsis <- function(ig, a=1, b=1.5) {
    N <- length(V(ig))
    # circle angles
    ang <- seq(0, 2*pi, length.out=(N+1))[-1]
    x <- a*cos(ang)
    y <- b*sin(ang)
    cbind(x,y)
}
plotdetailed <- function(phi, weights = NULL, main = "", stimuli = NULL, #layoutType = layout.graphopt,
		node.color="grey",node.size1=30, node.size2=7,
		edge.width=1,edge.arrowsize=0.5,
		layout=layout.circle, pdf=NULL, pointsize=12,
		edge.width.inhib=1.5, plot.legend=TRUE, label.cex=1, vlabel.cex=0.6, tk=FALSE,
		fontsize=20, rescale=TRUE)
{
	#stopifnot(require(igraph))
	edge.color="black"
	edge.lty="solid"
	phix <- phi
	phix[which(phix == 2)] <- 1
	ig <- graph.adjacency(phix)
	ig.nodes <- as.matrix(print.igraph.vs(V(ig))) # keggid
	if(tk) {
		tkplot(ig, vertex.label=V(ig)$name)
		return(ig)
	}
	vertex.color <- rep(node.color, length(ig.nodes))
	if(!is.null(stimuli)) {
		snodes <- unique(names(unlist(stimuli)))
		vertex.color[match(snodes, ig.nodes)] <- "red"
	}
	#browser()
	#inhib.coord <- t(sapply(which(phi==2),coord,mat=phi))
	inhib.coord <- which(phi==2,arr.ind=TRUE)
	if(!all(phix==0)) {
		ig.edges <- gsub(" ","",as.matrix(print.igraph.es(E(ig))))
		if(length(inhib.coord)!=0) {
			inhib.froms <- rownames(phi)[inhib.coord[,1]]
			inhib.tos <- rownames(phi)[inhib.coord[,2]]
			inhib.ed <- paste(inhib.froms,inhib.tos,sep="->")
			edge.color <- rep(edge.color, length(ig.edges))
			edge.lty <- rep(edge.lty, length(ig.edges))
			#edge.color[match(inhib.ed, ig.edges)] <- "darkgreen"
			edge.lty[match(inhib.ed, ig.edges)] <- "dashed"
			edge.width <- rep(edge.width, length(ig.edges))
			edge.width[match(inhib.ed, ig.edges)] <- edge.width.inhib
		}
		if(!is.null(weights)) {
			wind <- t(sapply(ig.edges,function(x) strsplit(x,split="->")[[1]]))
			edge.labels <- apply(wind, 1, function(x,ww) ww[x[1],x[2]], ww=weights)
			E(ig)$label <- edge.labels
			E(ig)$label.cex <- label.cex
			E(ig)$label.color <- "darkred"
		}
		## make reciprocated edges be curved
		ecurved <- rep(0,length(ig.edges))
		ed <- which(phix!=0,arr.ind=TRUE)
		erecipn <- NULL
		for(ei in 1:nrow(ed)) {
			edd <- ed[ei,]
			eddg <- phix[edd[1],edd[2]]
			eddgr <- phix[edd[2],edd[1]]
			if(eddgr!=0)
				erecipn <- c(erecipn,paste(rownames(phix)[edd[1]],colnames(phix)[edd[2]],sep="->"))
		}
		ecurved[match(erecipn,ig.edges)] <- 0.3
		E(ig)$curved <- ecurved
		
		## attributes
		print("Setting graph attributes...")
		E(ig)$color <- edge.color
		E(ig)$lty <- edge.lty
		E(ig)$width <- edge.width
		E(ig)$arrow.size <- edge.arrowsize
		#E(ig)$curved <- -0.3
		#E(ig)$curved <- -1 #TRUE
		V(ig)$label.cex <- vlabel.cex #0.6
	}
	V(ig)$color <- vertex.color
	V(ig)$shape <- rep("crectangle",length(ig.nodes))
	V(ig)$size <- rep(node.size1,length(ig.nodes))
	V(ig)$size2 <- rep(node.size2,length(ig.nodes))
	
	## set the layout		
	ig$layout <- layout
	
	print("..and plot.")
	if(!is.null(pdf))
		pdf(pdf,pointsize=pointsize)
	par(mar=c(1,1,1,1)) #, oma=c(0,0,0,0))
	plot(ig, vertex.label=V(ig)$name, main=main, rescale=rescale)
	#plot(ig, vertex.label=V(ig)$name, layout = layout, main=main)
	if(plot.legend) {
		legend("bottomright",legend=c("activation","inhibition"), col=c("black","black"), bty="n",
			lty=c("solid","dashed"), lwd=c(2,2), inset=c(0,-0.1), xpd=NA)
	}
	par(mar=c(5,4,4,2))
	if(!is.null(pdf))
		dev.off()
	invisible(list(ig=ig, layout=layout))
}

