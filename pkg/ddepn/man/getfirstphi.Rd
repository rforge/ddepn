\name{getfirstphi}
\Rdversion{1.1}
\alias{getfirstphi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Helper function for netga.
}
\description{
  Internal helper. Sample an initial network.
}
\usage{
getfirstphi(x, dat, stimuli, V, tps, reps, hmmiterations,
 lambda=NULL, B=NULL, Z=NULL, fanin=4, gam=NULL, it=NULL, 
 K=NULL, priortype="none", scale_lik=FALSE,
 allow.stim.off=TRUE, implementation="C")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A candidate network, given as adjacency matrix.}
  \item{dat}{The data matrix.}
  \item{stimuli}{The input stimuli list.}
  \item{V}{The names of the nodes.}
  \item{tps}{The time points.}
  \item{reps}{The number of replicates.}
  \item{hmmiterations}{The maximum number of iterations for the HMM.}
  \item{lambda}{The Prior influence strength in the laplace prior.}
  \item{B}{The Prior information matrix. Corresponds to prior edge probabilities in the final network.}
  \item{Z}{The normalisation factor for the prior distribution.}
  \item{fanin}{Integer: maximal indegree for nodes.}
  \item{gam}{Prior influence strength in scalefree prior.}
  \item{it}{Number of iterations to generate background distribution in scalefree prior.}
  \item{k}{Proportionality factor in scalefree prior}
  \item{priortype}{Character. Type of prior to be used.} 
  \item{scale_lik}{Boolean. Scale liklihood according to number of data points.}
  \item{allow.stim.off}{Boolean. Allow the stimuli to become inactive at some point. See also \code{\link{ddepn}}.}
  \item{implementation}{String. One of \code{"C","R","R_globalest","C_globalest"}. Different implementations of the HMM in 
    \code{\link{perform.hmmsearch}}. If \code{"R"}, the original pure R-implementation is 
    used, if \code{"C"}, a ported C-implementation is used. If \code{"R_globalest"}, an experimental version of
    the parameter estimation is used in the HMM, \code{"C_globalest"} is the C-port of 
    this version. See \code{details} for a description.}
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
\code{\link{ddepn}}
}
%\examples{
%## lala
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

