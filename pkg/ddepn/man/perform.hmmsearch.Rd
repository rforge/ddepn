\name{perform.hmmsearch}
\alias{perform.hmmsearch}
\alias{getE}
\alias{replicatecolumns}
\alias{perform.hmmsearch_R}
\alias{perform.hmmsearch_C}
\alias{perform.hmmsearch_globalest}
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
}
\details{
\code{perform.hmmsearch} is called whenever an HMM is to be run. Depending on the global variable \code{IMPLEMENTATION}, it decides 
upon the implementation of the HMM to use. \code{IMPLEMENTATION} can be set by passing command line argument \code{implementation} 
in the function call to \code{\link{ddepn}}. See the documentation there for a description.
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

