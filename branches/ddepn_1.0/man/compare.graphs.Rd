\name{compare.graphs.tc}
\Rdversion{1.1}
\alias{compare.graphs.tc}

\title{
compare.graphs.tc
}
\description{
Given an original and inferred network, count number of true/false positives and negatives. 
Calculate sensitivity and specificity measures.
}
\usage{
compare.graphs.tc(O,M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{O}{Original graph. Stored as adjacency matrix.}
  \item{M}{Inferred graph. Stored as adjaceny matrix.}
}
\details{
Counts the number of \eqn{tp} as: activations in \eqn{M} and \eqn{O} + inhibitions in \eqn{M} and \eqn{O}
Counts the number of \eqn{fp} as: activations or inhibitions in \eqn{M} and no edge in \eqn{O}
Counts the number of \eqn{tn} as: no edge in \eqn{M} and no edge in \eqn{O}
Counts the number of \eqn{fn} as: no edge in \eqn{M} and activation or inhibition in \eqn{O}
}
\value{
Vector of comparison measures:
\item{tp}{True positives}
\item{fp}{False positives}
\item{tn}{True negatives}
\item{fn}{False negatives}
\item{sn}{Sensitivity: \eqn{sn = tp/(tp+fn)}}
\item{sp}{Specificity: \eqn{sp = tn/(tn+fp)}}
\item{prec}{Precision: \eqn{prec = tp / (tp + fp)}}
\item{f1}{F1: \eqn{2*prec*sn / (prec + sn)}}
}
\author{
Christian Bender
}
\examples{
\dontrun{
  M <- matrix(sample(c(0,1,2),9,replace=TRUE),nrow=3,ncol=3)
  O <- matrix(sample(c(0,1,2),9,replace=TRUE),nrow=3,ncol=3)
  compare.graphs.tc(M,O)
}
}

%\keyword{ ~coordinates }
%\keyword{ ~matrix }% __ONLY ONE__ keyword per line
