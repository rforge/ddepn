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
getfirstphi(x, datx, stimuli, V, tps, reps, maxiter, lambda=NULL, B=NULL, Z=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	A candidate network, given as adjacency matrix.
  }
  \item{datx}{
	The data matrix.
  }
  \item{stimuli}{
	The input stimuli list.
  }
  \item{V}{
	The names of the nodes.
  }
  \item{tps}{
	The time points.
  }
  \item{reps}{
	The number of replicates.
  }
  \item{maxiter}{
	The maximum number of iterations for the HMM.
  }
  \item{lambda}{
  	The Prior influence strength.
  }
  \item{B}{
  	The Prior information matrix.
  }
  \item{Z}{
  	The normalisation factor for the prior distribution.
  }
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

