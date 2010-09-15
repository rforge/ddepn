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
posterior(phi, L, lambda=NULL, B=NULL, Z=NULL, gam=NULL, it=NULL, K=NULL)
pi(i,gam,N)
pgs(phi,gam,K=0.8,it=500)
}
\arguments{
  \item{phi}{	The candidate network.  }
  \item{L}{	The likelihood computed by \code{\link{likl}}.  }
  \item{lambda}{	Laplace prior hyperparameter describing the prior influence strength.  }
  \item{B}{	Laplace prior probability matrix.  }
  \item{Z}{	Laplace prior normalisation factor for the prior.  }
  \item{gam}{  Scale-free prior	degree distribution coefficient: \eqn{P(k) ~ k^gam}  }
  \item{N}{  	Number of nodes  }
  \item{K}{  	Scale-free prior scaling factor/Strength  }
  \item{it}{  Scale-free prior number of iterations for prior sampling.  }
}
\details{
Computes the posterior density from the likelihood and the prior information lambda, B and Z.
}
\value{
A double containing the posterior density.
}
\references{
\emph{Laplace prior}\cr
Froehlich et. al. 2007, Large scale statistical inference of signaling pathways from RNAi
and microarray data.\cr
Wehrli and Husmeier 2007,  Reconstructing gene regulatory networks with bayesian networks by combining
  expression data with multiple sources of prior knowledge\cr
  
\emph{Scale free prior}\cr
Kamimura and Shimodaira, A Scale-free Prior over Graph Structures for Bayesian Inference of Gene Networks
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

