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
		hmmiterations=5, phiasis=FALSE, lambda=NULL, B=NULL, Z=NULL, fanin=4,
		gam=NULL, it=NULL, K=NULL, priortype="none", scale_lik=FALSE,
		 allow.stim.off=TRUE, implementation="C")
initialphi(dat, phi, stimuli, Lmax, thetax, gammax, gammaposs,
		tps, reps, antibodies, n=100, multicores=FALSE, lambda=NULL, 
		B=NULL, Z=NULL, gam=NULL, it=NULL, K=NULL, priortype="none", 
		scale_lik=FALSE, allow.stim.off=TRUE, implementation="C")
}
\arguments{
  \item{dat}{The data matrix}
  \item{phi}{The initial network. Can be NULL.}
  \item{stimuli}{The stimuli list.}
  \item{hmmiterations}{Integer specifying the maximum number of iterations in the statespacesearch}
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
  \item{lambda}{The Prior influence strength in the laplace prior.}
  \item{B}{The Prior information matrix. Corresponds to prior edge probabilities in the final network.}
  \item{Z}{The normalisation factor for the prior distribution.}
  \item{fanin}{Integer: maximal indegree for nodes.}
  \item{gam}{Prior influence strength in scalefree prior.}
  \item{it}{Number of iterations to generate background distribution in scalefree prior.}
  \item{k}{Scaling factor in scalefree prior}
  \item{priortype}{Character. Type of prior to be used.}
  \item{scale_lik}{Boolean.Scale likelihood according to number of data points from which the 
    overall likelihood is calculated.}
  \item{allow.stim.off}{Boolean. Allow the stimuli to become inactive at some point. See also \code{\link{ddepn}}.}
  \item{implementation}{String. One of \code{"C","R","R_globalest","C_globalest"}. Different implementations of the HMM in 
    \code{\link{perform.hmmsearch}}. If \code{"R"}, the original pure R-implementation is 
    used, if \code{"C"}, a ported C-implementation is used. If \code{"R_globalest"}, an experimental version of
    the parameter estimation is used in the HMM, \code{"C_globalest"} is the C-port of 
    this version. See \code{details} for a description.}
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
