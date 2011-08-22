\name{plot_profiles}
\alias{plot_profiles}
\alias{get_theta_consensus}
\alias{get_gamma_consensus}
\alias{heatmapcolors}

\title{
Plot the data and inferred model parameters.
}
\description{
Creates a time profile for each experiment and protein in the data set. A spline fit is 
shown representing the trend of the intensities over time. Finally, model parameter 
estimations obtained in the inference are laid over the plots to show the active or 
passive states. 
}
\usage{
plot_profiles(ret, log=FALSE, ord=NULL, mfrow=c(4,4), 
              plotcurves=TRUE, plothist=TRUE, selection.criterion="aic")

get_theta_consensus(ret)

get_gamma_consensus(ret)

heatmapcolors(dat,ncol,lowcol="green",highcol="red",middlecol="white")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ret}{
    List or numeric data matrix. Either an object returned by \code{\link{ddepn}} with arguments \code{inference="netga"} 
    for Genetic Algorithm optimisation or \code{inference="mcmc"} for inhibMCMC sampling. Can also be a raw data matrix,
    in which case no parameter estimates are shown.
}
  \item{log}{
    Boolean. Transform data to log scale.
}
  \item{ord}{
    Vector of Strings. Optional, defines a node order for the plots.
}
  \item{mfrow}{
    Vector of numerics. Two values giving the number of rows and columns in which the output plot should be
    arranged.
}
\item{plotcurves}{
    Boolean. Show the data profiles and fitted splines or not.
}
  \item{plothist}{
    Boolean. Show additional histogram representation of the data or not.
}
  \item{selection.criterion}{
    For the spline fits to the data, choose a model selection criterion from \code{"aic"},
    \code{"bic"} or \code{"pvalue"} (for Akaike Information Criterion, Bayesian Information Criterion or 
    p-value selection. This will affect how many degrees of freedom are used for the spline fit. 
}
  \item{dat}{
    The data to be represented by the colour values. Is used to retrieve the range of values for which 
    the colour palette is created. 
}
  \item{ncol}{
    Integer. The number of colours in the panel.
}
  \item{lowcol}{
    A colour string, either keyword or hex-representation. Defining the colour used for the passive state.
}
  \item{middlecol}{
    A colour string, either keyword or hex-representation. Defining the colour in the middle of the palette.
}
  \item{highcol}{
    A colour string, either keyword or hex-representation. Defining the colour used for the active state.
}
}
\details{
Plots for each protein and experiment the data points along the time axis (as boxplots). Fits a smoothing spline to 
the data (see \code{\link{bestgam}}), using model selection criterion \code{selection.criterion} for deciding 
on the number of degrees of freedom. If \code{ret} is a list (from GA or inhibMCMC inference), the estimated
model parameters for the Gaussian active and passive distributions are extracted and shown in the timecourse 
plots (active: red, passive:green). If \code{ret} is a numeric matrix, only the data will be plotted. Boxplots
are coloured using a diverging colourpanel (created by \code{heatmapcolors} from green to red, indicating the state of the node at the 
respective time point.
}
\value{
none
}
\author{
Christian Bender
}
\note{
Needs package \emph{gam} for the spline fitting.
}
\seealso{
\code{\link{bestgam}}
}
\examples{
\dontrun{
## load package
library(ddepn)
library(multicore)
library(gam)

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
            
## perform netga inference, using 4 CPU cores
ret2 <- ddepn(dataset$datx, phiorig=phit, maxiterations=20, p=15, q=0.3, m=0.8, plotresults=FALSE, inference="netga", usebics=FALSE, priortype="laplaceinhib", lambda=0.01, B=B, multicores=TRUE, cores=4)

## plot the data and fits
## mcmc
plot_profiles(ret, mfrow=c(3,4))
## netga
plot_profiles(ret2, mfrow=c(3,4))
## data only
plot_profiles(ret2$dat, mfrow=c(3,4))



}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
