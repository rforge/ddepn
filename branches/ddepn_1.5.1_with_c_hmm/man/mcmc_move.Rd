\name{mcmc_move}
\alias{mcmc_move}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mcmc_move - Perform MCMC step
}
\description{
One step in an MCMC sampling. Changes exactly one edge in the actual network.
}
\usage{
mcmc_move(bestmodel, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bestmodel}{
	The current model.
  }
  \item{type}{
	String: The type of the move. One in addactivation, addinhibition, delete, reverse and switchtype. 
  }
}
\details{
TODO
}
\value{
TODO
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
\seealso{
\code{\link{mcmc_ddepn}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{ internal }

