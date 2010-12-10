\name{mcmc_performance}
\alias{mcmc_performance}
\alias{trapezoid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mcmc_performance - Performance plot for MCMC
}
\description{
Creates ROC and AUC for MCMC sampling.
}
\usage{
mcmc_performance(lst, plot = TRUE)
trapezoid(sp,sn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{
	The result list from an MCMC run.
  }
  \item{plot}{
  Boolean. If TRUE, show a plot of the ROC curve.
  }
}
\details{
\code{trapezoid} calculates the (approximate) area under curve (AUC) of two vectors of SN and SP values.
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

