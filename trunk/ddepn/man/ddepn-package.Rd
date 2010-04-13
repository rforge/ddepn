\name{ddepn-package}
\Rdversion{1.1}
\alias{ddepn-package}
\docType{package}
\title{
Dynamic Deterministic Effects Propagation Networks: Infer signalling networks from RPPA-Data.
}
\description{
Uses reverse phase protein array data to infer the signalling patterns among the measured proteins
under different stimuli/inhibitions. Each stimulus/inhibition is included as a separate node 
and edges to the measured proteins are inferred. Models the signal flow through the network and 
defines a set of system state transitions that can be reached in the given network. Searches for
the optimal path through these possible system states using an HMM. Estimation of Gaussian distributions
for active and passive states for each node depending on the activity state of the protein is done. Network search
is performed in a genetic algorithm.
}
\details{
\tabular{ll}{
Package: \tab ddepn\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-01-08\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Christian Bender

Maintainer: Christian Bender <c.bender@dkfz-heidelberg.de>
}
%\references{
%TODO
%}
\keyword{ package }
\seealso{
\code{\link{ddepn}}, \code{\link{netga}}
}
\examples{
\dontrun{
## load package
library(ddepn)

## sample a network
n <- 6
signet <- signalnetwork(n=n, nstim=2, cstim=0, prop.inh=0.2)
phit <- signet$phi
stimuli <- signet$stimuli

## sample data
dataset <- makedata(phit, stimuli, mu.bg=1200, sd.bg=400, mu.signal.a=2000, sd.signal.a=1000)

## uniform prior matrix
B <- matrix(0.5, nrow=n, ncol=n, dimnames=dimnames(phit))

## Genetic algorithm
ret <- ddepn(dataset$datx, phiorig=phit, phi=NULL, stimuli=NULL, th=0.5, inference="netga",
              pdf=NULL, multicores=FALSE, maxiterations=30, p=15, q=0.3, m=0.8, P=NULL,
              usebics=TRUE, cores=10, lambda=NULL, B=NULL, maxiter=100)	
x11()
plotdetailed(ret$phi)
}
}
