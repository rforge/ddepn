\name{ddepn}
\Rdversion{1.1}
\alias{ddepn}
\title{
ddepn
}
\description{
Main function for ddepn-modelling. Takes a data matrix of timecourse measurements as argument
and infers the most likely network structure underlying the data.
}
\usage{
ddepn(dat, phiorig=NULL, phi=NULL, stimuli=NULL, th=0.5, inference="netga", pdf=NULL,
      multicores=FALSE, maxiterations=1000, p=500, q=0.3, m=0.8, P=NULL,
      usebics=TRUE, cores=2, lambda=NULL, B=NULL, maxiter=100, fanin = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
		%%     ~~Describe \code{groupName} here~~
		Matrix of double values. The data matrix to be used. Contains antibody measurements in
		the rows and experiments (T timepoints in each R replicates) in the columns. Each
		experiment is labeled	by the respective perturbation in the column name. See section Details
		for an example.S}
	\item{phiorig}{
		Adjacency matrix. Original network used for comparison to the inferred net. Entries
		can be either 0, 1 or 2, for \emph{no edge}, \emph{activation} or \emph{inhibition},
		respectively. NULL if no reference network is given.}
	\item{phi}{
		Adjacency matrix. Network to start the search with. Entries
    	can be either 0, 1 or 2, for \emph{no edge}, \emph{activation} or \emph{inhibition},
    	respectively. NULL if no start network should be initialised automatically.}	
	\item{stimuli}{
		List containing the different perturbations in the system. See section Details for an
		example. }
	\item{th}{Threshold for inclusion of an edge in the final network. If an edge occures 
    	more than \eqn{th*p} times in all individuals, it is included in the 
    	resulting network.}
	\item{inference}{String giving the type of network search. 
		\code{netga} Uses a genetic algorithm for network inference.
		\code{mcmc} MCMC sampling for network inference.}
	\item{pdf}{String describing a path. If given, store output plots here}
	\item{multicores}{TRUE for using multiple cores and parallelize the repeated network
		reconstruction. FALSE for standard sequential handling of all repeats (needs
		R-package multicore).}
	\item{maxiterations}{Integer, Maximum number of generations in
		  \code{\link{netga}} or maximum number of iterations in \code{\link{mcmc_ddepn}}.}
	\item{p}{Integer, number of individuals in
		  the population in \code{\link{netga}}.}
	\item{q}{Double in [0;1], selection (1-q) and crossover (q) rate in \code{\link{netga}}.}
	\item{m}{Double in [0;1], mutation rate in \code{\link{netga}}.}
	\item{P}{List containing initial population of networks
		  for \code{\link{netga}}. Set to NULL if start population should be generated automatically.}
	\item{usebics}{Use BIC statistic for model selection. (only for \code{\link{netga}})}
	\item{cores}{Number of cores to use in case of \code{multicores}=TRUE}
	\item{maxiter}{Maximum number of iterations in the HMM search.}
	\item{lambda}{The Prior influence strength.}
  	\item{B}{The Prior information matrix.}
  	\item{fanin}{Integer: maximal indegree for nodes.}
}
\details{
  \describe{
    \item{dat}{
      Data matrix. Rows correspond to measured proteins/genes etc. Columns contain
      all experiments, i.e. separate perturbations. Each experiment consists of T
      timepoints and each timepoints was measured in R replicates. Example:
      \tabular{lcccccccc}{
         \tab EGF\_1\tab EGF\_1 \tab EGF\_T \tab EGF\_T \tab HRG\_1 \tab HRG\_1 \tab HRG\_T \tab HRG\_T\cr
         EGF \tab 0 \tab  0 \tab  0  \tab 0 \tab  0 \tab  0  \tab 0 \tab  0 \cr
         HRG \tab 0 \tab  0 \tab  0  \tab 0  \tab 0 \tab 0 \tab  0  \tab 0 \cr
         AKT \tab 1.45 \tab 1.8 \tab  0.99 \tab 1.6 \tab 1.78 \tab 1.8  \tab 1.56 \tab 1.58 \cr
         ERK \tab 1.33 \tab 1.7 \tab  1.57 \tab 1.3  \tab 0.68 \tab 0.34 \tab 0.62 \tab 0.47 \cr
         MEK \tab 0.45 \tab 0.8 \tab  0.99 \tab 0.6  \tab 0.78 \tab 0.8  \tab 0.56 \tab 0.58 \cr  
      }
    }
    \item{stimuli}{
      A list of stimuli describes the input perturbations, for which the activity states are
      fixed to be active. For the above example data matrix, stimuli would be
      \code{list(list(EGF=1),list(HRG=2))}. Also combinations of input stimuli are allowed, 
      e.g.  \code{list(list(EGF=1),list(HRG=2), list(EGF=1,HRG=2))} for a data set containing
      two distinct stimuli for EGF and HRG, and one stimulus EGF+HRG together.
    }
  }
}
\value{
	A list of the following elements:
	\item{phi}{The inferred network}
	\item{phi.activation.count}{Counts how often an edge is an activation in the population.}
	\item{phi.inhibition.count}{Counts how often an edge is an inhibition in the population.}
	\item{conf.act}{phi.activation.count/p}
	\item{conf.inh}{phi.inhibition.count/p}
	\item{weights}{Matrix: each entry is the maximum of the conf.act/conf.inh entries. I.e. this describes
		the support for an edge in the final network.}
	\item{weights.tc}{Similar to weights, but calculated ignoring the types of the edges.}
  \item{result}{A matrix containing result statistics for each network in the population:
          TP, FP, TN, FN, Sensitivity(SN), Specificity(SP), precision, F1. Only present if a 
          reference network phi.orig was provided in the function call.}	
  \item{P}{The population of networks that was inferred. NULL if inference="mcmc".}
	
	\item{phi.orig}{The reference network, if it was provided.}
	\item{L}{Likelihood}
	\item{aic}{Akaikes Information Criterion}
	\item{bic}{Bayesian Information Criterion}
	\item{dat}{The data}
	\item{theta}{The parameter matrix for the gaussians.}
	\item{gamma}{The state vectors for all timepoints and experiments}
	\item{tps}{Timepoints}
	\item{stimuli}{List of stimuli}
	\item{reps}{Number of replicates}
	\item{maxiter}{Maximum iterations to performed.}
}
\references{
%% ~put references to the literature/web site here ~
Froehlich et.al. 2009, Deterministic Effects Propagation Networks
  for Reconstructing Protein Signaling Networks from Multiple Interventions
Wehrli and Husmeier 2007,  Reconstructing gene regulatory networks with bayesian networks by combining
  expression data with multiple sources of prior knowledge
%Froehlich et.al. 2008, Estimating Large Scale Signaling Networks through Nested Effect Models with
%  Intervention Effects from Microarray Data
}
\author{
%%  ~~who you are~~
Christian Bender
}
\note{
TODO
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
TODO
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
              usebics=TRUE, cores=10, lambda=NULL, B=NULL, maxiter=100,fanin=2)	
## mcmc
ret <- ddepn(dataset$datx,phiorig=phit,phi=NULL,stimuli=stimuli,
      th=0.5, inference="mcmc", pdf=NULL, usebics=FALSE, multicores=FALSE,
      maxiterations=100, lambda=5, B=B, maxiter=50)

x11()
plotdetailed(ret$phi)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~network }
%\keyword{ ~greedy }% __ONLY ONE__ keyword per line
%\keyword{ ~simulatedannealing }
%\keyword{ ~Likelihood }
%\keyword{ ~statespacesearch }
