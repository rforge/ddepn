\name{plotdetailed}
\Rdversion{1.1}
\alias{plotdetailed}
\alias{get.labels}
\alias{get.arrowhead}
\alias{get.arrowtail}

\title{
plotdetailed
get.labels
get.arrowhead
get.arrowtail
}
\description{
\code{plotdetailed} uses Rgraphviz to plot a graph containing activations and inhibitions. Input is an
adjacency matrix with activation edge represented by a 1, inhibition by a 2.
}
\usage{
plotdetailed(phi,weights=NULL,main="",stimuli=NULL,
			layoutType = "dot", fontsize=20)
get.labels(phi)
get.arrowhead(phi)
get.arrowtail(phi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{The network to be plotted.}
  \item{weights}{Optional matrix of edge weights.}
  \item{main}{Optional character string containing main title for the plot.}
  \item{stimuli}{List of input stimuli. Are marked as red filled nodes in the plot.}
  \item{layoutType}{The graphviz layout for drawing the network.}
  \item{fontsize}{Fontsize of text in the plot.}
}
\details{
\code{get.labels}, \code{get.arrowhead} and \code{get.arrowtail} are helper functions for
the plot construction and usually not called directly.

The \code{stimuli} list describes the input stimuli or inhibitions. Each of these external treatment
nodes are drawn as red nodes in the graph. The format of the list should be:\\
Example: Assume two treatments \code{EGF} and \code{X}, a data matrix data:
      \tabular{lcccccccc}{
         \tab EGF_1\tab EGF_1 \tab EGF_2 \tab EGF_2 \tab EGF&X_1 \tab EGF&X_2 \tab EGF&X_2 \tab EGF&X_2\cr
         EGF \tab 0 \tab  0 \tab  0  \tab 0 \tab  0 \tab  0  \tab 0 \tab  0 \cr
         X   \tab 0 \tab  0 \tab  0  \tab 0  \tab 0 \tab 0 \tab  0  \tab 0 \cr
         AKT \tab 1.45 \tab 1.8 \tab  0.99 \tab 1.6 \tab 1.78 \tab 1.8  \tab 1.56 \tab 1.58 \cr  
      }
The stimuli list should be in the format:
\code{stimuli <- list(c(EGF=1),c(EGF=1,X=2))},
i.e. each element in the list corresponds to a treatment. Each treatment is a named numeric vector 
describing the row indices of the treatments in the data matrices. 
}
\value{
none
}
\author{
Christian Bender
}
\examples{
\dontrun{
library(ddepn)
mat <- matrix(sample(c(0,1,2),25,replace=TRUE), nrow=5,ncol=5, dimnames=list(LETTERS[1:5],LETTERS[1:5]))
weights <- matrix(sample(c(1:25),25,replace=TRUE), nrow=5,ncol=5, dimnames=list(LETTERS[1:5],LETTERS[1:5]))
stimuli <- list(list(A=1))
plotdetailed(mat,weights,"Some random graph",stimuli,layoutType="dot")
plotdetailed(mat,weights,"Some random graph",stimuli,layoutType="neato")
}
}
