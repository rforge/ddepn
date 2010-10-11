\name{zlambda}
\Rdversion{1.1}
\alias{zlambda}
\title{
Compute normalisation factor for prior information.
}
\description{
Compute normalisation factor for prior as defined in Wehrli/Husmeier 2007.
}
\usage{
zlambda(B, lambda)
}
\arguments{
  \item{B}{	The prior information matrix. }
  \item{lambda}{The strength of prior influence.  }
}
\details{
 Compute normalisation factor for Laplace prior as defined in Wehrli/Husmeier 2007.\cr 
 \eqn{Z = prod_v ( sum_pa(v) ( exp(-lambda * (sum_minpa(v)(1-B[v,m]) + sum_mnotinpa(v)(B[v,m])) ) ) )
   =                                           B[-v,pa(v)]}
}
\value{
Normalisation factor Z.
}
\references{
Wehrli and Husmeier 2007,  Reconstructing gene regulatory networks with bayesian networks by combining
  expression data with multiple sources of prior knowledge\cr
}
\author{
Christian Bender
}
\note{
Not used at the moment
}
\seealso{
TODO
}
\examples{
\dontrun{
library(ddepn)
B <- matrix(runif(16),nrow=4,ncol=4,dimnames=list(LETTERS[1:4],LETTERS[1:4]))
zlambda(B,2.2)
}
}
\keyword{ internal }
