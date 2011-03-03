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
  \item{gam}{ Scalefree prior degree distribution coefficient: \eqn{P(k) ~ k^gam} or
   exponent for difference term in laplaceinhib prior.  }
  \item{K}{   Scale-free prior scaling factor/Strength  }
  \item{it}{  Scale-free prior number of iterations for prior sampling.  }
  \item{priortype}{String. Either \code{uniform}, \code{laplaceinhib}, \code{laplace} or \code{scalefree}.
    \code{uniform} assumes uniform prior network distributions. 
    \code{laplaceinhib} calculates the difference between the candidate network and a 
    reference matrix containing edge probabilites, while both edge types (activation and 
    inhibition) are included. \code{laplace} is the same, except for ignoring the edge type. 
    \code{scalefree} calculates a probability of a network following a scale free network
    architecture. See the references for a detailed description of the priors.}
}
\details{
For the \emph{laplaceinhib} and \emph{laplace} prior types, the matrix \eqn{B} is of 
central importance. The matrix has the same dimensions as the network to be inferred,
each entry corresponding to a confidence in the existence of the respective edge. This
confidence can be aquired by using external pathway sources, e.g. the KEGG database. 
See the vignette for a description of how to get the prior matrices. No matter how
the confidence scores are obtained, there are two options, either use the \emph{laplaceinhib}
prior type, in which knowledge about the type of the edges is present in the 
external pathway source. Each confidence score for an edge, that is found as 
inhibiting edge in the reference pathways, is multiplied by \eqn{-1} to obtain a 
negative value for the inhibiting edges. The entries of the prior matrix \eqn{B} 
thus lie in the interval \eqn{[-1;1]}, where \eqn{-1} means strong confidence that 
an edge is an inhibition, \eqn{1} means strong confidence that the edge is an activation 
and \eqn{0} means that nothing is known about the presence or type of the edge.\cr

If no information on the type of the edges is available in the external data source,
\emph{priortype="laplace"} should be used, where the edge confidence ranges in the 
interval \eqn{[0;1]}, where \eqn{1} means strong confidence that the edge is present 
and \eqn{0} means that nothing is known about the presence of the edge.\cr  

Argument \code{gam} is used either as exponent in the scalefree prior, as it is described
in the reference, or in \code{laplaceinhib} and \code{laplace} as exponent in the following formula:\cr
\eqn{P(phi_{ij}|lambda,gam,B) = \frac{1}{2 \cdot lambda} exp(\frac{-|phi_{ij}-B_{ij}|^{gam}}{lambda})}\cr
It controls how strong the differences between an inferred edge and the probability
for seeing this edge in a reference set of networks are to be weighted. Defaults to 1,
if omitted. The prior curve rapidly decayes with increasing difference 
\eqn{|phi_{ij}-B_{ij}|}, while for \code{gam} larger than 1, the prior curve is
changed to an s-shaped curve with a plateau at the upper bound of
\eqn{P(phi_{ij}|lambda,gam,B)} and an exponential decay when a certain threshold of 
\eqn{|phi_{ij}-B_{ij}|} is reached.\cr

See also the references for a description of the priors.
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
