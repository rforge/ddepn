\name{mcmc_ddepn}
\alias{mcmc_ddepn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mcmc_ddepn - Perform MCMC sampling for DDEPN.
}
\description{
MCMC sampling for DDEPN. Takes an initial network and samples from the posterior.
}
\usage{
mcmc_ddepn(dat, phiorig = NULL, phi = NULL, stimuli = NULL, th = 0.8,
		 multicores = FALSE, pdf = NULL, maxiterations = 10000, 
		 usebics = FALSE, cores = 2, lambda = NULL, B = NULL, Z = NULL,
		 maxiter = 30, fanin = 4, gam = NULL, it = NULL, k = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{The data matrix.}
  \item{phiorig}{The reference network to compare to. Can be NULL.}
  \item{phi}{The start network. Empty if NULL.}
  \item{stimuli}{The stimuli list.}
  \item{th}{Threshold for inclusion of an edge in the final network.}
  \item{multicores}{Use multiple cores. Not used here.}
  \item{pdf}{File to which the network should be drawn.}
  \item{maxiterations}{Integer. Maximum number of MCMC iterations.}
  \item{usebics}{Use bics for model selection.}
  \item{cores}{Not used here.}
  \item{lambda}{Strength of prior influence.}
  \item{B}{Prior information matrix.}
  \item{Z}{Normalisation factor for prior.}
  \item{maxiter}{Maximum iterations in the HMM.}
  \item{fanin}{Integer: maximal indegree for nodes.}
  \item{gam}{Prior influence strength in sparsity prior.}
  \item{it}{Number of iterations to generate background distribution in sparsity prior.}
  \item{k}{Proportionality factor in sparsity prior}
}
\details{
MCMC is still under development. Use it with care.
}
\value{
  A list of the following elements:
  \item{phi}{The inferred network}
  \item{L}{Likelihood}
  \item{aic}{Akaikes Information Criterion}
  \item{bic}{Bayesian Information Criterion}
  \item{posterior}{Posterior probability}
  \item{dat}{The data}
  \item{theta}{The parameter matrix for the gaussians.}
  \item{gamma}{The state transition matrix.}
  \item{tps}{Timepoints}
  \item{stimuli}{List of stimuli}
  \item{reps}{Number of replicates}
  \item{maxiter}{Maximum iterations to performed.}
  \item{TSA}{deprecated}
  \item{Tt}{deprecated}  
  \item{lastmove}{Type of the last change that was performed.}
  \item{coords}{Position in the network where the last change was performed.}
  \item{lambda}{Prior strength parameter.}
  \item{B}{Prior matrix}
  \item{Z}{Prior normalisation factor.}
  \item{pegm}{Probability of performing the last move.}
  \item{pegmundo}{Probability of reverting the last move.}
  \item{nummoves}{Total number of possible moves in the current step}
  \item{stats}{Matrix of scores and statistics recorder during the MCMC run.}
  \item{freqa}{Counts how often an edge was an activation.}
  \item{freqi}{Counts how often an edge was an activation.}
}
\references{
TODO
}
\author{
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

## prior normalisation factor
lambda <- 5 
Z <- zlambda(B, lambda)

## network to start with
V <- rownames(dataset$datx)
phistart <- matrix(0, nrow=n, ncol=n, dimnames=list(V,V))

## now the sampling
ret <- mcmc_ddepn(dataset$datx, phiorig=phit, phi=phistart, stimuli=stimuli,
				th=0.5, multicores=FALSE, maxiterations=100,
				usebics=usebics, lambda=5, B=B, Z=Z, maxiter=30)
}
%ret <- ddepn(dataset$datx,phiorig=phit,phi=NULL,stimuli=stimuli,
%			th=th, inference="mcmc", pdf=NULL, usebics=FALSE,
%			maxiterations=maxiterations, lambda=lambda, B=B, maxiter=30)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
