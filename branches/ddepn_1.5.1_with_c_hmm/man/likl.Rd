\name{likl}
\Rdversion{1.1}
\alias{likl}
\title{
likl - Likelihood computation for given network, state-vectors and parameter matrix.
}
\description{
Given a data matrix \code{dat} and a matrix of state vectors \code{gammax} get the parameter
matrix \code{theta} and compute the Likelihood of the data \eqn{p(D|gammax,theta)} as sum of
all normal densities of each data point with the corresponding parameter chosen from theta. 
}
\usage{
likl(dat,gammax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data matrix as described in \link{ddepn}.}
  \item{gammax}{Matrix containing state vectors in the columns. Must have same dimensions as \code{dat}.}
}
\value{
	A list containing the likelihood matrix L and the parameter matrix theta.
	\item{L}{Matrix with dimensions dim(\code{dat}) containing the normal densities for
		 each data point given the state from gamma and the corresponding parameter from theta.
	}
	\item{thetax}{Parameter matrix containing the means and standard deviations for the nodes.
  		A node corresponds to a row in the data matrix. Each row contains a number of time points,
  		which are separated in active and passive states. For each state, one mean/sd parameter
  		pair is found in \code{theta}.
	}
}
\references{
TODO
}
\author{
Christian Bender
}
\seealso{
\code{\link{ddepn}}
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}
\keyword{ internal }

