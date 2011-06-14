\name{plot_mcmctraces}
\alias{plot_mcmctraces}

\title{
Trace plot for inhibMCMC posterior traces.
}
\description{
Plots the traces of posterior probabilities during an inhibMCMC run.
}
\usage{
plot_mcmctraces(ret, thin = 1)
}
\arguments{
  \item{ret}{
    List. Object returned by \code{\link{ddepn}} call, using argument 
    \code{inference="mcmc"} for inhibMCMC sampling.
}
  \item{thin}{
    The thinning interval between consecutive observations.
}
}
\details{
Plot the MCMC posterior traces of all runs during inhibMCMC inference. Uses a thinning interval 
defined in \code{thin}. 
}
\value{
none
}
\author{
Christian Bender
}
\note{
Uses the \code{coda} package and its mcmc object and plot functions.
}
\seealso{
\code{\link{ddepn}}
}
\examples{
\dontrun{
## load package
library(ddepn)
library(multicore)
library(coda)

## sample a network and data
set.seed(1234)
n <- 6
signet <- signalnetwork(n=n, nstim=2, cstim=0, prop.inh=0.2)
phit <- signet$phi
stimuli <- signet$stimuli
dataset <- makedata(phit, stimuli, mu.bg=1200, sd.bg=400, mu.signal.a=2000, sd.signal.a=1000)

## use original network as prior matrix
## reset all entries for inhibiting edges 
## to -1
B <- phit
B[B==2] <- -1

## perform inhibMCMC inference, using 4 CPU cores to get 4 MCMC chains 
ret <- ddepn(dataset$datx, phiorig=phit, maxiterations=300, burnin=50,
            plotresults=FALSE, inference="mcmc",
            usebics=FALSE, priortype="laplaceinhib", lambda=0.01, B=B,
            multicores=TRUE, cores=4)

plot_mcmctraces(ret,thin=1)
plot_mcmctraces(ret,thin=10)

}
}
