\name{kegggraphs}
\Rdversion{1.1}
\alias{kegggraphs}
\docType{data}
\title{
kegggraphs
}
\description{
Contains the KEGG Hierarchy. Downloaded around August 2009.
}
%\usage{data(kegggraphs)}
\format{List of 195, each Element is one KEGG pathway.}
\details{
The KEGG Hierarchy was downloaded and transformed into graphNEL objects
by using KEGGgraph and graph.
}
\source{
http://www.genome.jp/kegg/
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Kanehisa: Kyoto Encyclopaedia of Genes and Genomes
}
\examples{
\dontrun{
library(ddepn)
data(kegggraphs)
names(kegggraphs)
}
}
\keyword{datasets}
