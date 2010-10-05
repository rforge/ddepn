\name{ddepn-package}
\Rdversion{1.1}
\alias{ddepn-package}
\docType{package}
\title{
Dynamic Deterministic Effects Propagation Networks: Infer signalling networks from high througput array data.
}
\description{
Uses high throughput array data (e.g. from Reverse Phase Protein arrays or mRNA microarrays)
to infer regulatory or signalling relationships among the measured proteins after different stimuli/inhibitions.
Each stimulus/inhibition is included as a separate node and edges to the measured proteins are
inferred. The signal flow through the network is modelled in a boolean fashion and a set of reachable system states
is created. Afterwards, the optimal path through the possible system states is searched for using an HMM.
Estimation of Gaussian distributions for active and passive states for each node depending on the activity state of
the protein is done and a likelihood based score is calculated for a network. Network structure search is performed
in either a genetic algorithm that optimises a population of candidate networks or via a markov chain monte carlo sampling
that summarises the distribution over the network structures via sampling.
}
\details{
\tabular{ll}{
Package: \tab ddepn\cr
Type: \tab Package\cr
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

## Genetic algorithm, using BIC score as optimisation criterion
ret1 <- ddepn(dataset$datx, phiorig=phit,inference="netga",
              maxiterations=30, p=15, q=0.3, m=0.8,
              usebics=TRUE)	
x11()
plotdetailed(ret1$phi,stimuli=stimuli,fontsize=20)

## Genetic algorithm, using a laplaceinhib prior and the posterior probabilities as optimisation criterion
B <- matrix(0.5,nrow=n,ncol=n,dimnames=dimnames(phit))
ret2 <- ddepn(dataset$datx, phiorig=phit, inference="netga",
              maxiterations=30, p=15, q=0.3, m=0.8,
              usebics=FALSE, lambda=5, B=B, priortype="laplaceinhib") 
x11()
plotdetailed(ret2$phi,stimuli=ret2$stimuli,fontsize=20)

## MCMC sampling using a laplaceinhib prior
ret3 <- ddepn(dataset$datx,phiorig=phit, inference="mcmc",
              maxiterations=300, burnin=100,
              usebics=FALSE, lambda=2, B=B, priortype="laplaceinhib")

x11()
plotdetailed(ret3$samplings[[1]]$phi,stimuli=ret3$samplings[[1]]$stimuli)

## MCMC sampling using a scale free prior
ret4 <- ddepn(dataset$datx,phiorig=phit, inference="mcmc",
              maxiterations=300, burnin=100,
              usebics=FALSE, gam = 2.2, it = 500, K = 0.8, priortype="scalefree")
x11()
plotdetailed(ret4$samplings[[1]]$phi,stimuli=ret4$samplings[[1]]$stimuli)


}
}
