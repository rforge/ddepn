\name{getfirstphi}
\Rdversion{1.1}
\alias{getfirstphi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Helper function for netga.
}
\description{
  Internal helper for \code{\link{netga}}
}
\usage{
getfirstphi(x, dat, stimuli, V, tps, reps, maxiter,
 lambda=NULL, B=NULL, Z=NULL, fanin=4, gam=NULL, it=NULL, k=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A candidate network, given as adjacency matrix.}
  \item{dat}{The data matrix.}
  \item{stimuli}{The input stimuli list.}
  \item{V}{The names of the nodes.}
  \item{tps}{The time points.}
  \item{reps}{The number of replicates.}
  \item{maxiter}{The maximum number of iterations for the HMM.}
  \item{lambda}{The Prior influence strength in the laplace prior.}
  \item{B}{The Prior information matrix. Corresponds to prior edge probabilities in the final network.}
  \item{Z}{The normalisation factor for the prior distribution.}
  \item{fanin}{Integer: maximal indegree for nodes.}
  \item{gam}{Prior influence strength in sparsity prior.}
  \item{it}{Number of iterations to generate background distribution in sparsity prior.}
  \item{k}{Proportionality factor in sparsity prior}
}
\details{
Usually called internally.
}
\value{
TODO
}
\references{
TODO
}
\author{
Christian Bender
}
\note{
TODO
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
TODO
}
\examples{
## lala
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

