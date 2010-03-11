\name{perform.hmmsearch}
\alias{perform.hmmsearch}
\alias{getE}
\alias{replicatecolumns}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform.hmmsearch - Do system state optimisation. Usually not called directly.
}
\description{
Use an HMM to perform search for optimal system state sequence. \code{getE} and
\code{replicatecolumns} are two helper functions for computing the emission probabilities
and doing some internal matrix operations. \code{perform.hmmsearch} is usually called internally. 
}
\usage{
perform.hmmsearch(phi.n, bestmodel)
getE(x, datx, thetaprime)
replicatecolumns(mat, replicates=4)
}
%- maybe also 'usage' for other objects documented here.
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
TODO
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

