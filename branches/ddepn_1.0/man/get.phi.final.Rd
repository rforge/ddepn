\name{get.phi.final}
\alias{get.phi.final}
\title{
Construct final network from GA result.
}
\description{
Takes the output of \code{\link{netga}} and constructs a final network from the population.
}
\usage{
get.phi.final(lst, th = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{Output list from \code{\link{netga}}
  }
  \item{th}{Double in [0;1]. Threshold for inclusion of an edge into the final network.
  }
}
\details{
Takes the population P from the GA resultlist and returns the list with the element \emph{lst\$phi} 
replaced by the new final network. \emph{lst\$weights} only contains the weights
with \eqn{lst\$weights > th}.
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
## uniform prior matrix
B <- matrix(0.5, nrow=n, ncol=n, dimnames=dimnames(phit))
## Genetic algorithm
ret <- ddepn(dataset$datx, phiorig=phit, phi=NULL, stimuli=NULL, th=0.01, inference="netga",
              pdf=NULL, multicores=TRUE, maxiterations=30, p=15, q=0.3, m=0.8, P=NULL,
              usebics=TRUE, cores=10, lambda=NULL, B=NULL, maxiter=100) 

plotrepresult(ret,pdf=NULL)
ret2 <- get.phi.final(ret, th=0.9)
plotrepresult(ret2,pdf=NULL)
}
}

