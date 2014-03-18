\name{plot_overlap}
\alias{plot_overlap}
\title{
   Plot a reconstructed network in comparison with a prior network.
}
\description{
   Display a reconstructed network together with the provided prior network in the same plot.
}
\usage{
   plot_overlap(adjcore, adjX, overlap, symbols.plot, ig.layout, main, nameA, nameB, 
                node.size1=25, node.size2=7, stimnodes=NULL, customtitle=NULL,
                col_equal="black", col_new="blue", col_missing="red", vlabel.cex=1)
}
\arguments{
  \item{adjcore}{
		Adjacency matrix of the prior network (row target activates column target, 0: no edge, 1: activation, 2: inhibition).
 }
  \item{adjX}{
    Adjacency matrix of the reconstructed consensus network (row target activates column target, 0: no edge, 1: activation, 2: inhibition).
}
  \item{overlap}{
      Targets of prior and reconstructed network, i.e. row/column names of the adjacency matrices.
}
  \item{symbols.plot}{
    Vertex label, i.e. row/column names of the adjacency matrices.
}
  \item{ig.layout}{
  Layout matrix defining the positions of vertices in the plot. The row names are the targets and the column names are "x"" and "y"" as coordinates. The entries are numeric values.
}
  \item{main}{
  Plot subtitle.
 }
  \item{nameA}{
	Provide a name for the prior network, e.g. "Prior".
}
  \item{nameB}{
  Provide a name for the reconstructed network, e.g. "Reconstruction".
}
  \item{node.size1}{
  Node size parameter.
}
  \item{node.size2}{
  Node size parameter.
}
  \item{stimnodes}{
  Vector of stimuli names.
}
  \item{customtitle}{
  Plot title.
}
  \item{col_equal}{
  Colour for common edge in prior and reconstructed network.
}  
\item{col_new}{
 Colour for newly inferred edge in reconstructed network compared to prior.
}
 \item{col_missing}{
 Colour for deleted edge in reconstructed network compared to prior.
}
 \item{vlabel.cex}{
  Label size parameter.
}
}
\value{
	Plot of a reconstructed network in comparison with a prior network.
}

\author{
	Christian Bender
}

\examples{

	library(ddepn)
  ## simulate a network with 6 nodes and 2 distinct input stimuli
  set.seed(12345)
  signet <- signalnetwork(n=6, nstim=2, cstim=0, prop.inh=0.2)
  ## prior net
  prior <- signet$phi
  prior[5,3] <-  prior[2,1] <- 1
  prior[6,4] <- 0
  ## layout
  layoutMatrix<-matrix(0,nrow=6,ncol=2)
  rownames(layoutMatrix)<-rownames(prior)
  colnames(layoutMatrix)<-c("x","y")
  layoutMatrix[,"x"]<-rep(c(100,200),3)
  layoutMatrix[,"y"]<-c(100,100,200,200,300,300)
	plot_overlap(adjcore=prior, adjX=signet$phi, ig.layout=layoutMatrix,overlap=rownames(prior), nameA="Prior", nameB="Reconstruction", main="Test", symbols.plot=rownames(prior), stimnodes=names(unlist(signet$stimuli)))
}
