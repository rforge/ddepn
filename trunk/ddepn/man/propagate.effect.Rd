\name{propagate.effect.set}
\Rdversion{1.1}
\alias{propagate.effect.set}
\alias{propagate.effect.simple}

\title{
propagate.effect.set
propagate.effect.simple
}
\description{
Effect propagation through signalling network. Usually called internally.
}
\usage{
propagate.effect.set(phi,tps,stimuli,reps=1)
propagate.effect.simple(phi,tps,stimulus,reps=1)
}
\arguments{
  \item{phi}{Network}
  \item{tps}{Timepoints}
  \item{stimuli}{The input stimuli list.}
  \item{stimulus}{One input stimulus.}
  \item{reps}{Replicates}
}
\details{
\code{propagate.effect.set} performs effect propagation for multiple experiments, identified by
	the different stimuli. It calls \code{propagate.effect} separately for each stimulus. 
}
\value{
The state transition matrix containing all possible state vectors.
}
\author{
Christian Bender
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
\keyword{ internal }
