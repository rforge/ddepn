\name{posterior}
\alias{posterior}
\alias{pri}
\alias{pgs}
\title{
Posterior probability calculation.
}
\description{
Posterior probability calculation. Usually called internally.
}
\usage{
posterior(phi, L, lambda=NULL, B=NULL, Z=NULL, gam=NULL, it=NULL, 
		K=NULL, priortype="laplaceinhib")
pri(i,gam,N)
pgs(phi,gam,K=0.8,it=500)
}
\arguments{
  \item{phi}{	The candidate network.  }
  \item{L}{	The likelihood computed by \code{\link{likl}}.  }
  \item{lambda}{	Laplace prior hyperparameter describing the prior influence strength.  }
  \item{B}{	Laplace prior probability matrix.  }
  \item{Z}{	Laplace prior normalisation factor for the prior. (Not used at the moment.)  }
  \item{gam}{  Scale-free prior	degree distribution coefficient: \eqn{P(k) ~ k^gam}  }
  \item{N}{  	Number of nodes  }
  \item{K}{  	Scale-free prior scaling factor/Strength  }
  \item{it}{  Scale-free prior number of iterations for prior sampling.  }
  \item{priortype}{Character. One of \code{uniform}, \code{laplaceinhib}, \code{laplace} or
    \code{scalefree} for use of the respective prior type.} 
}
\details{
Computes the posterior density depending on \code{priortype}: \code{uniform} uses uniform prior, \code{laplaceinhib} and \code{laplace} use 
prior parameters lambda, gam, B and Z, and \code{scalefree} uses gam and K as prior parameters. See \code{\link{prior}} for a description of 
the prior models.
}
\value{
A double containing the posterior density.
}
\references{
\emph{Laplace prior}\cr
Froehlich et. al. 2007, Large scale statistical inference of signaling pathways from RNAi
and microarray data.
%\cr
%Wehrli and Husmeier 2007,  Reconstructing gene regulatory networks with bayesian networks by combining
%  expression data with multiple sources of prior knowledge\cr
  
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
## TODO
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

