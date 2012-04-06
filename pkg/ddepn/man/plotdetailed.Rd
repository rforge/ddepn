\name{plotdetailed}
\alias{plotdetailed}
\alias{layout.ellipsis}

\title{
plotdetailed, layout.ellipsis
}
\description{
\code{plotdetailed} uses igraph to plot a graph containing activations and inhibitions. Input is an
adjacency matrix with activation edge represented by a 1, inhibition by a 2.
}
\usage{
plotdetailed(phi, weights = NULL, main = "", stimuli = NULL,
  node.color = "grey", node.size1 = 30, node.size2 = 7, edge.width = 1, edge.arrowsize = 0.5,
  layout = layout.circle, pdf = NULL, pointsize = 12, edge.width.inhib = 1.5, plot.legend = TRUE,
  label.cex = 1, vlabel.cex = 0.6, tk = FALSE, fontsize=20, rescale=TRUE)

layout.ellipsis(ig, a=1, b=1.5)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{Adjacency matrix. The network to be plotted.}
  \item{weights}{Optional matrix of edge weights.}
  \item{main}{Optional character string containing main title for the plot.}
  \item{stimuli}{List of input stimuli. Are marked as red filled nodes in the plot.}
  %-\item{layoutType}{The graphviz layout for drawing the network.}
  \item{fontsize}{Deprecated. Still there for compatibility issues with \code{\link{plotdetailed_Rgraphviz}}. Use \code{cex.label}
    to control edge label size.}
  \item{node.color}{Define coloring of nodes.}
  \item{node.size1}{Width of a rectangular node.}
  \item{node.size2}{Height of a rectangular node.}
  \item{edge.width}{Width of edges.}
  \item{edge.arrowsize}{Size of edge ends.}
  \item{layout}{Graph layout algorithm or matrix with node coordinates.}
  \item{pdf}{Pdf output file name.}
  \item{pointsize}{Pointsize argument used in pdf device.}
  \item{edge.width.inhib}{Width of inhibitory edges.}
  \item{plot.legend}{Define if legend should be drawn or not.}
  \item{label.cex}{Edge label size factor.}
  \item{vlabel.cex}{Node label size factor.}
  \item{tk}{Use tk-plotting engine or not.}
  \item{rescale}{Rescaling coordinates or not.}
  \item{ig}{Igraph object.}
  \item{a}{Semi major axis of ellipse.}
  \item{b}{Semi minor axis of ellipse.}
}
\details{
Plot a graph encoded in adjacency matrix \code{phi} using the \code{igraph} library. 
Activation edges are drawn as solid lines, inhibition edges as dashed lines. Stimuli nodes can be 
colored red using the \code{stimuli} argument.

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
describing the row indices of the treatments in the data matrices. \cr


\code{layout.ellipsis} is another graph layout making elliptical graphs possible. Can be used as layout argument for 
  plotting. Use \code{rescale=FALSE} if using layout.ellipsis.
}
\value{
A list object holding the igraph object \code{ig} and the layout coordinates \code{layout}:
	\item{ig}{Igraph object created from the adjacency matrix.}
	\item{layout}{Matrix with two columns, holding the x and y coordinates of the elements in the graph.}
}
\references{
none
}
\author{
Christian Bender
}

\seealso{
\code{\link{plotdetailed_Rgraphviz}}
}
\examples{
\dontrun{
library(ddepn)
mat <- matrix(sample(c(0,1,2),25,replace=TRUE), nrow=5,ncol=5, dimnames=list(LETTERS[1:5],LETTERS[1:5]))
diag(mat) <- 0
weights <- matrix(sample(c(1:25),25,replace=TRUE), nrow=5,ncol=5, dimnames=list(LETTERS[1:5],LETTERS[1:5]))
stimuli <- list(list(A=1))
plotdetailed(mat,weights,"Some random graph",stimuli,layout=layout.graphopt)
layoutmat <- cbind(c(1,2,3,4,5), c(5,4,3,2,1))
plotdetailed(mat,weights,"Some random graph",stimuli,layout=layoutmat)
}
}
