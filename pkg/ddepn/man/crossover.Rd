\name{crossover}
\alias{crossover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ crossover - perform crossing over of two networks }
\description{
Perform a crossover of two networks at a random node and
 exchange all connected edges.
}
\usage{
crossover(phi1, phi2, fanin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi1}{ Adjacency matrix, the first network. }
  \item{phi2}{ Adjacency matrix, the second network. }
  \item{fanin}{ Max indegree of node. }
}
\details{
Chooses a random node and performs crossover of network \code{phi1}
 and \code{phi2}. Used for \code{\link{netga}}.
}
\value{
  A list containing the matrices after the crossover.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Christian Bender }
\note{ ~~further notes~~ 
}
\seealso{ \code{\link{netga}} }
%\examples{
%TODO
%}
\keyword{ internal }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
