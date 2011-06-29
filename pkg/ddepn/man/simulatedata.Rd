\name{simulatedata}
\Rdversion{1.1}
\alias{simulatedata}


\title{
simulatedata
}
\description{
Generate artificial timecourse data for network inference.
}
\usage{
simulatedata(phi, mu.bg=0, sd.bg=0.1,
		mu.signal.a=1, sd.signal.a=0.5,
		mu.signal.i=-1, sd.signal.i=0.5,
		stimulus=sample(nrow(phi),2),TT=10,R.t=4,R.b=3,
		plot=FALSE,stimuli=NULL, allow.stim.off=TRUE) 
}
\arguments{
  \item{phi}{The network for which data should be simulated.}
  \item{mu.bg}{mean for passive state.}
  \item{sd.bg}{sd for passive state.}
  \item{mu.signal.a}{mean for active state of type activation.}
  \item{sd.signal.a}{sd for active state of type activation.}
  \item{mu.signal.i}{mean for active state of type inhibition.}
  \item{sd.signal.i}{sd for active state of type inhibition.}
  \item{stimulus}{Where the network gets stimulated. Are set to 1 for the effect propagation.}
  \item{TT}{Number of timepoints.}
  \item{R.t}{Number of technical replicates.}
  \item{R.b}{Number of biological replicates.}
  \item{plot}{Should a plot be generated after data generation.}
  \item{stimuli}{List of input stimuli.}
  \item{allow.stim.off}{Boolean. Allow the stimuli to become inactive at some point. See also \code{\link{ddepn}}.}
}
\details{
TODO
}
\value{
Artificial datasets for a given network.
}
\author{
Christian Bender
}

\examples{
\dontrun{
library(ddepn)
n <- 8
phi <- matrix(sample(c(0,1,2),n*n,replace=TRUE),nrow=n,dimnames=list(LETTERS[1:n],LETTERS[1:n]))
simulatedata(phi, mu.bg=0, sd.bg=0.1,
		mu.signal.a=1, sd.signal.a=0.5,
		mu.signal.i=-1, sd.signal.i=0.5,
		stimulus=sample(nrow(phi),2),TT=10,R.t=4,R.b=3,
		plot=TRUE)
}
}

%\keyword{ ~state transition matrix }
%\keyword{ ~expand }% __ONLY ONE__ keyword per line
