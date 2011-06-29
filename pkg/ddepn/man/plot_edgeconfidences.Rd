\name{plot_edgeconfidences}
\alias{plot_edgeconfidences}
\alias{make_edge_df}

\title{
Boxplots of edge confidences obtained in inhibMCMC.
}
\description{
Plot the confidences for each edge, obtained from \emph{L} independent MCMC chains.
}
\usage{
plot_edgeconfidences(ret, start = 1, stop = NULL, act = "conf.act", inh = "conf.inh", ...)

make_edge_df(samp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ret}{
    List. Either the returned object of the \code{\link{ddepn}} call using the inhibMCMC inference
    (argument \code{inference="mcmc"}), or the sublist \code{samplings} from this object. 
}
  \item{start}{
    Numeric. Optional, defines the start of the subset of the iterations in the MCMC chains which are
    used to generate the boxplots.
}
  \item{stop}{
    Numeric. Optional, defines the end of the subset of the iterations in the MCMC chains which are
    used to generate the boxplots.
}
  \item{act}{
    String. Defines the statistic to be used from the MCMC results. One of \code{"conf.act"} or \code{"freqa"},
    for activation confidence or frequency of activations.
}
  \item{inh}{
    String. Defines the statistic to be used from the MCMC results. One of \code{"conf.inh"} or \code{"freqi"},
    for inhibition confidence or frequency of inhibitions.
}
  \item{samp}{
    List. Contains the sampling runs from each MCMC chain. Note that the number of chains \eqn{L:=length(samp)}.
}
  \item{cex.axis}{
    Numeric. Scale factor for the axis labels.
}
  \item{\dots}{
    Further plotting options passed to the lattice \code{panel.bwplot} function.
}

}
\details{
Create a summary plot of edge confidences or counts over \emph{L} MCMC runs. Assume that \emph{N} is the number 
of nodes in the inferred network. Each panel in the summary plot contains \emph{N} boxes showing the activation 
and \emph{N} boxes showing the inhibition confidences/frequencies. Activation boxes are shown in blue, inhibition 
boxes in red. The column name in each panel defines the source node, from which an edge originates. The panel name 
denotes the destination node to which the edge points.
}
\value{
none
}
\author{
Christian Bender
}
\examples{
\dontrun{
## load package
library(ddepn)
library(multicore)

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

plot_edgeconfidences(ret, act="conf.act", inh="conf.inh",pch="|")

}
}
