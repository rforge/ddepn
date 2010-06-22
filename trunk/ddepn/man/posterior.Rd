\name{posterior}
\alias{posterior}
\alias{pi}
\alias{pgs}
\title{
Posterior probability calculation.
}
\description{
Posterior probability calculation. Usually called internally.
}
\usage{
posterior(phi, L, lambda, B, Z, gam=2.2, it=500, K=0.8)
pi(i,gam,N)
pgs(phi,gam,K=0.8,it=500)
}
\arguments{
  \item{phi}{
	The candidate network.
  }
  \item{L}{
	The likelihood computed by \code{\link{likl}}.
  }
  \item{lambda}{
	Weight for the prior influence.
  }
  \item{B}{
	The prior probability matrix.
  }
  \item{Z}{
	The normalisation factor for the prior.
  }
  \item{gam}{
  	degree distribution coefficient: \eqn{P(k) ~ k^gam}
  }
  \item{N}{
  	Number of nodes
  }
  \item{K}{
  	Scaling factor/Strength
  }
  \item{it}{
  	Number of iterations for prior sampling.
  }
}
\details{
Computes the posterior density from the likelihood and the prior information lambda, B and Z.
}
\value{
A double containing the posterior density.
}
\references{
Froehlich et.al. 2007, Wehrli & Husmeier 2007
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
## lala
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

