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
		plot=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{The network for which data should be simulated}
  \item{mu.bg}{mean for passive state}
  \item{sd.bg}{sd for passive state}
  \item{mu.signal.a}{mean for active state of type activation}
  \item{sd.signal.a}{sd for active state of type activation}
  \item{mu.signal.i}{mean for active state of type inhibition}
  \item{sd.signal.i}{sd for active state of type inhibition}
  %\item{mu.bg, sd.bg, mu.signal.a, sd.signal.a, mu.signal.i, sd.signal.i}{The parameters for 
  %	the normal distributions. Given for the background or 'passive' and the 'active' states,
  %	i.e. for activation/inhibition states.}
  \item{stimulus}{Where the network gets stimulated. Are set to 1 for the effect propagation.}
  %\item{stimuli}{List of stimuli.}
  \item{TT}{Number of timepoints}
  \item{R.t}{Number of technical replicates}
  \item{R.b}{Number of biological replicates}
  \item{plot}{Should a plot be generated after data generation}
  %\item{gammax}{State transition matrix}
  %\item{npairs}{Number of gamma/data pairs to be generated}
  %\item{partitions}{A vector containing a shift pattern for the gamma states in each row. Can 
  %	be given explicitely or is calculated internally (maybe very slow depending on the size
  %	of the network}
  %\item{pdf}{Where the output pictures should be saved}
  %\item{n}{number of nodes in the network}
  %\item{detailed}{TRUE for edges being modeled as activation/inhibition separately, FALSE for 
  %	only one edge type}
  %\item{nodes}{Character vector containing the node names}
  %\item{kegggraphslist}{A list containing the KEGG database}
  %\item{minEdges}{minimum number of edges in the network}
  %\item{must.contain.inhibitions}{TRUE if the network must contain at least one inhibition}
  %\item{onlysignaling}{Boolean, if TRUE, only signaling pathways and additionally Cell Cycle, 
  %	Apoptosis and p53 Signaling are used for subgraph sampling}
  %\item{datx}{The data matrix}
  %\item{kg}{A graphNEL object to sample nodes from}
  %\item{knodes}{Character vector, a subset of nodenames from kg}
  %\item{...}{Further arguments to be passed}
  %\item{mat}{The matrix to be plotted}
  %\item{name}{A plot title}
}
\details{
TODO
%\code{get.data} is an internal function and should not be called directly.
%\code{plot.matrix} is an internal function and should not be called directly.
%\code{get.partitions} is an internal function and should not be called directly.
}
\value{
Artificial datasets for a given network.
}
\author{
Christian Bender
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

phi <- matrix(sample(c(0,1,2),9,replace=TRUE),nrow=3,dimnames=list(LETTERS[1:3],LETTERS[1:3]))
simulatedata(phi, mu.bg=0, sd.bg=0.1,
		mu.signal.a=1, sd.signal.a=0.5,
		mu.signal.i=-1, sd.signal.i=0.5,
		stimulus=sample(nrow(phi),2),TT=10,R.t=4,R.b=3,
		plot=FALSE)
}

%\keyword{ ~state transition matrix }
%\keyword{ ~expand }% __ONLY ONE__ keyword per line
