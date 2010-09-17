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
compare.graphs.tc(phiorig,phi,ignore.type=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phiorig}{Original graph. Stored as adjacency matrix.}
  \item{phi}{Inferred graph. Stored as adjaceny matrix.}
  \item{ignore.type}{Boolean. If TRUE, only compare edge abundance, if FALSE, also take into account
    the type of the edge for counting numbers of true and false edges.}
}
\details{
Counts the number of \eqn{tp} as: activations in \eqn{phi} and \eqn{phiorig} + inhibitions in \eqn{phi} and \eqn{phiorig}
Counts the number of \eqn{fp} as: activations or inhibitions in \eqn{phi} and no edge in \eqn{phiorig}
Counts the number of \eqn{tn} as: no edge in \eqn{phi} and no edge in \eqn{phiorig}
Counts the number of \eqn{fn} as: no edge in \eqn{phi} and activation or inhibition in \eqn{phiorig}
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
  library(ddepn)
  phi <- matrix(sample(c(0,1,2),9,replace=TRUE),nrow=3,ncol=3)
  phiorig <- matrix(sample(c(0,1,2),9,replace=TRUE),nrow=3,ncol=3)
  compare.graphs.tc(phi,phiorig)
}
}

%\keyword{ ~coordinates }
%\keyword{ ~matrix }% __ONLY ONE__ keyword per line
