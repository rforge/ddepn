\name{get.phi.final}
\alias{get.phi.final}
\alias{get.phi.final.mcmc}
\title{
Construct final network from GA or mcmc result.
}
\description{
\code{get.phi.final} takes the output of \code{\link{netga}} and constructs a final network
from the population, \code{get.phi.final.mcmc} takes a list containing the returned lists of 
\code{\link{mcmc_ddepn}} and calculates a final network for each result list.
}
\usage{
get.phi.final(lst, th = 0.8)
get.phi.final.mcmc(retlist,maxiterations,prob=.333,qu=.99999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{Output list from \code{\link{netga}}  }
  \item{th}{Double in [0;1]. Threshold for inclusion of an edge into the final network.  }
  \item{retlist}{A list returned by MCMC inference. }
  \item{maxiterations}{Integer. Number of MCMC sampling iterations. }
  \item{prob}{Double \eqn{\in [0;1]}.Success probability for binomial density.}
  \item{qu}{Double \eqn{\in [0;1]}. Quantile of the binomial distribution, used as significance cutoff
    for edge inclusion.}
}
\details{
\code{get.phi.final}
Takes the population P from the GA resultlist and returns the list with the element \emph{lst$phi} 
replaced by the new final network. \emph{lst$weights} only contains the weights
with \eqn{lst$weights > th}.
\code{get.phi.final.mcmc}
Takes a list of MCMC samplings returned by \code{\link{mcmc_ddepn}} and extracts a final network for
each result. It is assumed that any type of edge is expected to be found with probability 1/3=0.33, thus
for \code{maxiterations} samplings we expect the probability of seeing exactly \code{k} edges of a certain
type follows a binomial distribution \code{B(maxiterations, prob, k)}. An edge is included, if the 
probability of seeing more than the observed number of occurrences is less or equal than \code{qu},
i.e. \eqn{P(k>K_{observed}) \leq qu}.
}
\value{
Result list as in \code{\link{netga}}, with replaced \emph{phi} and \emph{weights} elements.
}
\author{
Christian Bender
}

\seealso{
\code{\link{netga}}
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
## Genetic algorithm
ret1 <- ddepn(dataset$datx, phiorig=phit, inference="netga",
              maxiterations=30, p=15, q=0.3, m=0.8, P=NULL,
              usebics=TRUE) 

plotrepresult(ret1,outfile=NULL)
ret2 <- get.phi.final(ret1, th=0.9)
plotrepresult(ret2,outfile=NULL)

## mcmc
maxiterations <- 300
## use original network as prior matrix
## reset all entries for inhibiting edges 
## to -1
B <- phit
B[B==2] <- -1

ret3 <- ddepn(dataset$datx,phiorig=phit,
      inference="mcmc", usebics=FALSE,
      maxiterations=maxiterations, burnin=100, lambda=0.01, B=B)

plotrepresult(ret3$samplings[[1]],outfile=NULL)
ret4 <- get.phi.final(ret3$samplings[[1]],th=0.9)
plotrepresult(ret4,outfile=NULL)

}
}

