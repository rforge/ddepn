\name{gelman_diag}
\alias{gelman_diag}

\title{
Convergence check by Gelman's potential scale reduction.
}
\description{
Compute Gelman's potential scale reduction to estimate convergence in inhibMCMC.
}
\usage{
gelman_diag(ret)
}

\arguments{
  \item{ret}{List. Object returned by \code{\link{ddepn}}, when inhibMCMC was used as 
    inference method.
}
}
\details{
The function is a wrapper calling the \code{gelman.diag} function from package \code{coda}.
}
\value{
The result of the \code{gelman.diag} function.
}
\references{
Gelman, A., Carlin, JB., Stern, HS., Rubin, DB.: Bayesian Data Analysis, 2nd edition, 
  Chapman & Hall/CRC, chapter 11.6, pp294-295.
}
\note{
Note that multiple chains have to be performed to use this diagnostic.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
coda package.
}
\examples{
\dontrun{
library(ddepn)
set.seed(12345)
n <- 6
signet <- signalnetwork(n=n, nstim=2, cstim=0, prop.inh=0.2)
net <- signet$phi
stimuli <- signet$stimuli
weights <- signet$weights

## sample data
dataset <- makedata(net, stimuli, mu.bg=1200, sd.bg=400, mu.signal.a=2000, sd.signal.a=1000)
data <- dataset$datx

# use the original network as prior probability matrix
B <- net
B[B==2] <- -1
# construct a prior matrix with uniform probabilities for each edge
if(require(multicore)) {
  ret <- ddepn(data, phiorig=net, inference="mcmc",
    maxiterations=3000, burnin=1000,
    usebics=FALSE, lambda=0.01, B=B,
    multicores=TRUE, cores=4, priortype="laplaceinhib")
}

## now produce the convergence diagnostic
gelman_diag(ret)


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
