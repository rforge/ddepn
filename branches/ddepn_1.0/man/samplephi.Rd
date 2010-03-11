\name{samplephi}
\Rdversion{1.1}
\alias{samplephi}
\alias{initialphi}

\title{
sample random initial networks
}
\description{
Internal function.
}
\usage{
samplephi(phi,stimuli, antibodies, tps, reps, dat, searchstatespace=FALSE,
		maxiter=5, phiasis=FALSE, lambda=NULL, B=B, Z=Z)
initialphi(dat, phi, stimuli, Lmax, thetax, gammax, gammaposs,
		tps, reps, antibodies, n=100, multicores=FALSE, lambda=NULL, B=NULL, Z=NULL)
}
\arguments{
  \item{dat}{The data matrix}
  \item{phi}{The initial network. Can be NULL.}
  \item{stimuli}{The stimuli list.}
  \item{maxiter}{Integer specifying the maximum number of iterations in the statespacesearch}
  \item{multicores}{TRUE for using multiple cores.}
  \item{Lmax}{Likelihood score}
  \item{thetax}{Parameter matrix}
  \item{gammax}{Statespacematrix}
  \item{gammaposs}{All possible state vectors}
  \item{tps}{The timepoints}
  \item{reps}{Number of replicates}
  \item{antibodies}{Character vector of protein names in the network.}
  \item{searchstatespace}{Do statespacesearch in initial sampling?}
  \item{phiasis}{TRUE for taking a given phi matrix as sample}
}
\details{
Used for initializing random networks. Called internally.
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

\seealso{
TODO
}
\examples{
##TODO
}
\keyword{ internal }
