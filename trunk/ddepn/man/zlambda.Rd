\name{zlambda}
\alias{zlambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute normalisation factor for prior information.
}
\description{
Compute normalisation factor for sparsity prior as defined in Wehrli/Husmeier 2007.
}
\usage{
zlambda(B, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{
  	The prior information matrix.
  }
  \item{lambda}{
	The strength of prior influence.
  }
}
\details{
 Compute normalisation factor for sparsity prior as defined in Wehrli/Husmeier 2007.
 models prior knowledge as network energy that includes a boundary for maximal incoming 
 edges (fan.in).  
 Z = prod_v ( sum_pa(v) ( exp(-lambda * (sum_minpa(v)(1-B[v,m]) + sum_mnotinpa(v)(B[v,m])) ) ) )
   =                                           B[-v,pa(v)]
}
\value{
Normalisation factor Z.
}
\references{
Wehrli, Husmeier 2007
}
\author{
Christian Bender
}
\note{
TODO
}
\seealso{
TODO
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
\keyword{ internal }
