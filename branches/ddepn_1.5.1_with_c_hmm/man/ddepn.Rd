\name{ddepn}
\Rdversion{1.1}
\alias{ddepn}
\title{
ddepn
}
\description{
Main function for DDEPN modelling. Takes a data matrix containing
 longitudinal measurements as argument and infers a network
 structure underlying the data using either a genetic algorithm
 or MCMC sampling.
}
\usage{
ddepn(dat, phiorig=NULL, phi=NULL, th=0.5, inference="netga",
      outfile=NULL, multicores=FALSE, maxiterations=1000,
      p=500, q=0.3, m=0.8, P=NULL,
      usebics=TRUE, cores=2, priortype="laplaceinhib",
      lambda=NULL, B=NULL, samplelambda=TRUE,
      hmmiterations=100, fanin=4,
      gam=NULL,it=NULL,K=NULL,quantL=.5,quantBIC=.5,
      debug=FALSE, burnin=1000, thin=FALSE)
}
\arguments{
  \item{dat}{
		Matrix of double values. The data matrix to be used. Contains
		antibody measurements in the rows and experiments (T timepoints
		in each R replicates) in the columns. Each experiment is
		labeled	by the respective perturbation in the column name.
		See section Details for an example.}
	\item{phiorig}{
		Adjacency matrix. Reference network used for comparison to the
		inferred net. Entries can be either 0, 1 or 2, for \emph{no edge},
		\emph{activation} or \emph{inhibition}, respectively. NULL if no
		reference network is given.}
	\item{phi}{
		Adjacency matrix. Seed network to start the search. Entries
    	can be either 0, 1 or 2, for \emph{no edge}, \emph{activation}
    	or \emph{inhibition}, respectively. NULL if no start network
    	should be given, but initialised automatically.}	
	\item{th}{Threshold for inclusion of an edge in the final network
	 (for \code{netga}). If an edge occurs more than \eqn{th*p} times
	 in all individuals, it is included in the resulting network.}
	\item{inference}{String. Giving the type of network search.\cr 
		\code{netga} Uses a genetic algorithm for network inference.
		\code{mcmc} MCMC sampling for network inference.}
	\item{outfile}{String. Output path for plotting. NULL if plotting
	 should be done to the display.}
	\item{multicores}{Boolean. TRUE for using multiple cores and
	 parallelise the network reconstruction. In case of \code{netga}
	 the HMMs for each individual in the population are distributed
	 on multiple cores. In case of \code{mcmc}, several independent
	 MCMC runs are started, each on a separate core. FALSE for
	 standard calculation on only one core (needs R-package multicore).}
	\item{maxiterations}{Integer, Maximum number of generations in
		  \code{\link{netga}} or maximum number of iterations in
		  \code{\link{mcmc_ddepn}}.}
	\item{p}{Integer, number of individuals in the population in
	 \code{\link{netga}}.}
	\item{q}{Double \eqn{\in [0;1]}, selection (1-q) and crossover (q)
	 rate in \code{\link{netga}}.}
	\item{m}{Double \eqn{\in [0;1]}, mutation rate in \code{\link{netga}}.}
	\item{P}{List containing an initial population of networks
		  for \code{\link{netga}}. Set to NULL if start population should
		  be generated automatically.}
	\item{usebics}{Use BIC statistic for model selection (only for
	 \code{\link{netga}}).}
	\item{cores}{Number of cores to use in case of \code{multicores}=TRUE.
	 For \code{netga}, the parallel calculations of the HMMs are
	 distributed on \code{cores} cores, for \code{mcmc} \code{cores} 
	 independent MCMC runs are started.}
	\item{hmmiterations}{Integer. Maximum number of iterations in the HMM
	 search.}
	\item{lambda}{The Prior influence hyperparameter for the laplace prior.}
  \item{B}{The Prior information matrix. See \code{\link{prior}} for details.}
  \item{fanin}{Integer: maximal indegree for each node.}
  \item{gam}{Prior influence strength for scalefree prior. Also used as exponent
    in \code{laplaceinhib} prior: see \code{\link{prior}} for details. }
  \item{it}{Number of iterations to generate the background distribution
   for scalefree prior.}
  \item{K}{Proportionality factor for scalefree prior.}
   \item{quantL}{Quantile of Population Likelihood/Posterior, used as
    selection threshold in \code{netga}. Note that the Likelihood or
    Posterior have to be maximised, so all networks with a
    likelihood/posterior \emph{greater} than this threshold are selected.}
  \item{quantBIC}{Quantile of Population BIC, used as selection
   threshold in \code{netga}. Note that the BIC is minimised, so all
    networks with BIC \emph{less} than the threshold are selected.}
  \item{samplelambda}{Boolean. TRUE if laplace prior hyperparameter
   \code{lambda} should be sampled, too. FALSE if fixed \code{lambda}
   is to be used. Unused if \code{lambda} or \code{B} equal zero.}
  \item{debug}{Boolean. If TRUE, print extra informational messages.}
  \item{burnin}{Integer. Specifies the number of iterations used as
   burnin phase for \code{\link{mcmc_ddepn}}.}
  \item{priortype}{Character. One of \code{none}, \code{uniform}, \code{laplaceinhib},
   \code{laplace} or \code{scalefree} for use of the respective prior
   type. Ignored if \code{usebics=TRUE} for \code{netga}. For \code{netga},
   \code{usebics=FALSE, priortype="none"} means optimising the likelihood 
   directly. This is equivalent to setting \code{usebics=FALSE, priortype="uniform"}.
   For \code{mcmc_ddepn}, \code{priortype="none"} is not allowed. Use
   \code{priortype="uniform"}instead. 
   \code{laplaceinhib} uses prior information for edges with two types
   (activation/inhibition), \code{laplace} ignores the edge type. Useful
   if only knowledge about the presence of an edge is available, but not 
   about its type. \code{scalefree} assumes scale-free network architectures.} 
  \item{thin}{Boolean. If TRUE, make sure that the MCMC return objects are
   shortened to at most 10000 iterations. Defaults to FALSE.}
}
\details{
  \describe{
    \item{dat}{
      Data matrix. Rows correspond to measured proteins/genes etc.
      Columns contain all experiments, i.e. separate perturbations.
      Each experiment consists of T time points and each time point is
      assumed to be measured in R replicates. The time is indicated as a
      numeric value, separated by an underscore in the column name.
      Example:
      \tabular{lcccccccc}{
      \tab EGF_1\tab EGF_1 \tab EGF_2 \tab EGF_2 \tab EGF&X_1 \tab EGF&X_2 \tab EGF&X_2 \tab EGF&X_2\cr
      EGF \tab 0 \tab  0 \tab  0  \tab 0 \tab  0 \tab  0  \tab 0 \tab  0 \cr
      X   \tab 0 \tab  0 \tab  0  \tab 0  \tab 0 \tab 0 \tab  0  \tab 0 \cr
      AKT \tab 1.45 \tab 1.8 \tab  0.99 \tab 1.6 \tab 1.78 \tab 1.8  \tab 1.56 \tab 1.58 \cr
      ERK \tab 1.33 \tab 1.7 \tab  1.57 \tab 1.3  \tab 0.68 \tab 0.34 \tab 0.62 \tab 0.47 \cr
      MEK \tab 0.45 \tab 0.8 \tab  0.99 \tab 0.6  \tab 0.78 \tab 0.8  \tab 0.56 \tab 0.58 \cr  
      }
    }
  }
}
\value{
	\emph{For \code{netga}, a list containing the following elements:}
	\item{dat}{Double matrix. The data matrix.}
	\item{phi.activation.count}{Integer. Counts how often an edge is an
	 activation in the population.}
	\item{phi.inhibition.count}{Integer. Counts how often an edge is an
	 inhibition in the population.}
	\item{phi.orig}{Adjacency matrix. The reference network, if it was
	 provided.}
	\item{phi}{Adjacency matrix. The inferred network}
	\item{weights}{Matrix. Each entry is the maximum of the conf.act/conf.inh
	 entries. I.e. this describes
    the support for an edge in the final network.}
  \item{weights.tc}{Matrix. Similar to weights, but calculated ignoring the
   types of the edges.}
  \item{stats}{Matrix. Contains result statistics for each network in the
   population: TP, FP, TN, FN, Sensitivity(SN), Specificity(SP), precision, F1.
   Only present if a reference network \code{phi.orig} was provided in the
   function call to \code{ddepn}.} 
	\item{conf.act}{Matrix. Calculated as phi.activation.count/p}
	\item{conf.inh}{Matrix. Calculated as phi.inhibition.count/p}
	\item{stimuli}{List. The given list of input stimuli.}
  \item{P}{List. The population of networks that was inferred, i.e. the
   return list of \code{\link{netga}}.}
  \item{scorestats}{Matrix. Contains traces of the scores during the genetic
  algorithm. See \code{\link{netga}}.}
	\cr
	\cr
	\emph{For \code{mcmc}, a list containing two elements:}
	\item{samplings}{List. Contains all sampling runs. Each sampling run itself
	 is a list as obtained via \code{\link{mcmc_ddepn}}.}
	\item{ltraces}{Matrix. Contains the posterior traces, each trace stored in
	 one column of the matrix.} 
}
\references{
\emph{DDEPN}\cr
Bender et. al. 2010: Dynamic deterministic effects propagation networks: learning 
 signalling pathways from longitudinal protein array data; Bioinformatics,
 Vol. 26(18), pp. i596-i602\cr
 
%% ~put references to the literature/web site here ~
%Froehlich et.al. 2009, Deterministic Effects Propagation Networks
%  for Reconstructing Protein Signaling Networks from Multiple Interventions\cr 

\emph{Laplace prior}\cr
Froehlich et. al. 2007, Large scale statistical inference of signaling pathways from RNAi
and microarray data; BMC Bioinformatics, Vol. 8(11), pp. 386ff\cr

%\cr
%Werhli and Husmeier 2007,  Reconstructing gene regulatory networks with bayesian networks by combining
%  expression data with multiple sources of prior knowledge\cr
  
\emph{Scale free prior}\cr
Kamimura and Shimodaira, A Scale-free Prior over Graph Structures for Bayesian
 Inference of Gene Networks

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
dataset <- makedata(phit, stimuli, mu.bg=1200, sd.bg=400,
                    mu.signal.a=2000, sd.signal.a=1000)

## use original network as prior matrix
## reset all entries for inhibiting edges 
## to -1
B <- phit
B[B==2] <- -1

## Genetic algorithm, no prior
ret1 <- ddepn(dataset$datx, phiorig=phit, inference="netga",
              maxiterations=30, p=15, q=0.3, m=0.8,
              usebics=TRUE)	
x11()
plotdetailed(ret1$phi,stimuli=ret1$stimuli)
              
## mcmc, laplaceinhib prior
ret2 <- ddepn(dataset$datx,phiorig=phit, inference="mcmc",
              maxiterations=300, burnin=100,
              usebics=FALSE, lambda=0.01, B=B, priortype="laplaceinhib") 
      
x11()
plotdetailed(ret2$samplings[[1]]$phi,stimuli=ret2$samplings[[1]]$stimuli)

## use mcmc with multiple cores, i.e. perform two independent runs
## requires package multicore and, of course multiple cores in the hardware
## use the original net as prior
 if(require(multicore)) {
 	ret3 <- ddepn(dataset$datx,phiorig=phit, inference="mcmc",
                multicores=TRUE, cores=2,
                maxiterations=300, burnin=100,
                usebics=FALSE, lambda=0.01, B=B, priortype="laplaceinhib")
 }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~network }
%\keyword{ ~greedy }% __ONLY ONE__ keyword per line
%\keyword{ ~simulatedannealing }
%\keyword{ ~Likelihood }
%\keyword{ ~statespacesearch }
