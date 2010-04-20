\name{propagate.effect.set}
\Rdversion{1.1}
\alias{propagate.effect.set}
\alias{propagate.effect.simple}
\alias{bin2dec}

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
bin2dec(x)
}
\arguments{
  \item{phi}{Network}
  \item{tps}{Timepoints}
  \item{stimuli}{The input stimuli list.}
  \item{stimulus}{One input stimulus.}
  \item{reps}{Replicates}
  \item{x}{binary vector}
}
\details{
\code{propagate.effect.set} performs effect propagation for multiple experiments, identified by
	the different stimuli. It calls \code{propagate.effect.simple} separately for each stimulus.
\code{bin2dec} takes a 0/1 vector and converts it to a decimal number. 
}
\value{
\code{propagate.effect.set} returns the state transition matrix containing all
possible state vectors.
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
