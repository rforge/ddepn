\name{coord}
\Rdversion{1.1}
\alias{coord}

\title{
coord - compute coordinates for a given index.
}
\description{
Given an index k and a matrix mat, compute the coordinate pair (x,y) at which k is found in mat.
}
\usage{
coord(k,mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{Integer value.}
  \item{mat}{Matrix.}
}
\details{
TODO
}
\value{
Integer pair (x,y)
}
\author{
Christian Bender
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
mat <- matrix(0, nrow=5,ncol=10)
k <- 14
coord(k,mat)
}
}

\keyword{ internal }
%\keyword{ ~matrix }% __ONLY ONE__ keyword per line
