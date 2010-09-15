\name{signalnetwork}
\Rdversion{1.1}
\alias{signalnetwork}
\title{
Sample a random signalling network.
}
\description{
For \eqn{n} nodes, sample a network with \eqn{nstim} stimuli and \eqn{cstim} combinatorial stimuli.
}
\usage{
signalnetwork(n=10, nstim=2, cstim=0, prop.inh=.2,plot=F,gamma=1,B=NULL)
}
\arguments{
  \item{n}{Integer. Number of nodes.}
  \item{nstim}{Integer. Number of stimuli.}
  \item{cstim}{Integer. Number of combinatorial stimuli.}
  \item{prop.inh}{Proportion (in [0;1]) of the number of activating edges to be
    included as inhibiting edges in the network.}
  \item{plot}{Boolean. If TRUE, a plot of the generated graph is drawn.}
  \item{gamma}{Double. Strength of power law decay. Used for simulating the
    number of outgoing edges.}
  \item{B}{The prior edge probability matrix}
}
\details{
	Simulates an artificial signalling network. Starts at \eqn{nstim} random stimuli and 
	selects random children, to which activation edges are drawn. These children are the new
	stimuli and the procedure is repeated until all nodes were reached by activating edges. 
	Finally, \eqn{prop.inh*numedges} inhibiting edges are added randomly. The number of
	stimuli combinations \eqn{cstim} is limited by \eqn{sum_{k=2}^n {k \choose n}}.
	If defined, B gives a matrix containing prior probabilities for each possible edge in 
	the network.
}
\value{
	List containing the adjacency list \eqn{phi} and the list of all stimuli.
}
\author{
Christian Bender
}
\seealso{
\code{\link{simulatedata}}
}
\examples{
\dontrun{
signalnetwork(n=10, nstim=4, cstim=4, prop.inh=.4, plot=TRUE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
