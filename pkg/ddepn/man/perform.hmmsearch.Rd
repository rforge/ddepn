\name{perform.hmmsearch}
\alias{perform.hmmsearch}
\alias{getE}
\alias{replicatecolumns}
\alias{perform.hmmsearch_R}
\alias{perform.hmmsearch_C}
\alias{perform.hmmsearch_globalest}
\alias{perform.hmmsearch_C_globalest}
\alias{updateA}
\alias{viterbi}
\alias{is_consistent}
\alias{checkC}
\title{
perform.hmmsearch - Do system state optimisation using a Hidden Markov Model. Usually not called directly.
}
\description{
Use an HMM to perform search for optimal system state sequence. \code{getE} and
\code{replicatecolumns} are two helper functions for computing the emission probabilities
and doing some internal matrix operations. \code{perform.hmmsearch} is usually called internally. 
}
\usage{
perform.hmmsearch(phi.n, bestmodel)

perform.hmmsearch_R(phi.n, bestmodel)

perform.hmmsearch_C(phi.n, bestmodel)

perform.hmmsearch_globalest(phi.n, bestmodel)

getE(x, datx, thetaprime)

replicatecolumns(mat, replicates=4)

is_consistent(phi, gammax)

checkC(phi, gammax, ch, pa)

updateA(maxima.ind, pseudocount, M, A)

viterbi(gammaposs, datx, thetaprime, T, R, M, A)

}
\arguments{
  \item{phi.n}{Adjacency matrix describing the candidate network.
  }
  \item{datx}{The data matrix.
  }
  \item{bestmodel}{Model object containing the current model.
  }
  \item{x}{State vector.
  }
  \item{thetaprime}{Model parameters for gaussians.
  }
  \item{mat}{Matrix object.
  }
  \item{replicates}{How often should each column be replicated.
  }
  \item{phi}{An adjacency matrix.}
  \item{gammax}{A state transition matrix.}
  \item{ch}{Integer denoting a child node.}
  \item{pa}{Integer denoting a parent node.}
  \item{maxima.ind}{Vector of integer. Series of state indices describing the 
    optimal path through the system states.}
  \item{pseudocount}{Integer. The pseudocount value used in the viterbi update of the state 
    transition matrix.}
  \item{M}{Integer. Number of states in A. Used for normalisation of the transition probability update.}
  \item{A}{The transition probability matrix.}
  \item{gammaposs}{The system state matrix containing \code{M} system states.}
  \item{T}{Integer. The number of time points.}
  \item{R}{Integer. The number of replicates.}
}
\details{
\code{perform.hmmsearch} is called whenever an HMM is to be run. Depending on the global variable \code{IMPLEMENTATION}, it decides 
upon the implementation of the HMM to use. \code{IMPLEMENTATION} can be set by passing command line argument \code{implementation} 
in the function call to \code{\link{ddepn}}. See the documentation there for a description.\cr

\code{is_consistent} checks whether a state transition matrix is consistent with the structure. This means, for 
instance, that a child node \code{ch} can never be activated whenever a parent \code{pa} is active and 
connected via an inhibition edge. The function returns the number of found inconsitencies and some extra 
information. It is called in implementation \code{R_globalest} to improve the assignment of active and
passive states. 
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{ internal }

