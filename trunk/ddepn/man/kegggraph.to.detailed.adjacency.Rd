\name{kegggraph.to.detailed.adjacency}
\Rdversion{1.1}
\alias{kegggraph.to.detailed.adjacency}

\title{
kegggraph.to.detailed.adjacency
}
\description{
Converts a graphNEL object extracted by KEGGgraph to an adjacency
matrix containing inhibitions and activations.
}
\usage{
kegggraph.to.detailed.adjacency(gR)
}
\arguments{
  \item{gR}{A KEGGgraph graphviz object.}
}
\details{
Can be used to convert the objects in \code{\link{kegggraphs}} to 
adjacency matrices.
}
\value{
An adjacency matrix.
}
\author{
Christian Bender
}

\examples{
\dontrun{
  data(kegggraphs)
  adj <- kegggraph.to.detailed.adjacency(kegggraphs[[1]]$g)
}
}

%\keyword{ internal }
%\keyword{ ~matrix }% __ONLY ONE__ keyword per line
