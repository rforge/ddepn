\name{prior}
\alias{prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate a structure prior. Usually called internally.
}
\description{
Provides three types of structure priors. \code{laplaceinhib} and \code{laplace} penalise the 
difference between the actual network and a reference network. \code{scalefree} penalises high
node degrees.
}
\usage{
prior(phi, lambda = NULL, B = NULL, Z = NULL,
   gam = NULL, it = NULL, K = NULL, priortype = "laplaceinhib")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{ The candidate network.  }
  \item{lambda}{  Laplace prior hyperparameter describing the prior influence strength.  }
  \item{B}{ Laplace prior probability matrix.  }
  \item{Z}{ Laplace prior normalisation factor for the prior. (Not used at the moment.)  }
  \item{gam}{ Scalefree prior degree distribution coefficient: \eqn{P(k) ~ k^gam} or exponent for difference term in 
  	laplaceinhib prior.  }
  \item{K}{   Scale-free prior scaling factor/Strength  }
  \item{it}{  Scale-free prior number of iterations for prior sampling.  }
  \item{priortype}{String. Either \code{laplaceinhib}, \code{laplace} or \code{scalefree}. 
    \code{laplaceinhib} calculates the difference between the candidate network and a 
    reference matrix containing edge probabilites, while both edge types (activation and 
    inhibition) are included. \code{laplace} is the same, except for ignoring the edge type. 
    \code{scalefree} calculates a probability of a network following a scale free network
    architecture. See the references for a detailed description of the priors.}
}
\details{
A description of the priors used is given in the references.\cr
The normal use of \code{gam} is for the scalefree prior, as it is described in the reference. If \code{laplaceinhib}
is used, \code{gam} is used as exponent in the following formula:\cr
\eqn{P(phi_{ij}|lambda,gam,B) = \frac{1}{2 \cdot lambda} exp(\frac{-|phi_{ij}-B_{ij}|^{gam}}{lambda})}\cr
It controls how strong the differences between an inferred edge and the probability for seeing this edge in a reference set
of networks are to be weighted. Defaults to 2, if omitted, which leads to weak penalisation of small differences between 
inferred and reference interaction, and strong penalisation for big differences.
%\section{laplaceinhib}{TODO}
%\section{laplace}{TODO}
%\section{scalefree}{TODO}
}
\value{
Returns a double for the prior probability of network structure \code{phi}.
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
\code{\link{ddepn}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
