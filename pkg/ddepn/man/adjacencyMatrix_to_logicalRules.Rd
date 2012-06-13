\name{adjacencyMatrix_to_logicalRules}
\alias{adjacencyMatrix_to_logicalRules}
\title{
Convert adjacency matrix to logical rules.
}
\description{
Convert adjacency matrices of \code{\link{ddepn}} consensus nets to logical activation/inhibition rules of network component, i.e. target, wiring for further analysis with the \cite{BoolNet-package}, e.g. for perturbation simulations, according to its \cite{loadNetwork} function format.}
\usage{
adjacencyMatrix_to_logicalRules(adjMatrix, outfile)
}
\arguments{
  \item{adjMatrix}{Adjacency matrix as resulting from \code{\link{create_signetwork}} for example.}
  \item{outfile}{Path to txt file where to store the output.}
}
\details{
The input adjacency matrix must have the network components, e.g. the protein names, as row names which must be identical to the column names. As typical for the consensus net obtained from \code{\link{ddepn}} network reconstructions, the following number coding holds for the matrix:
\itemize{
\item 1: target in row activates target in column,
\item 2: target in row inhibits target in column,
\item 0: no relation of net components.
}
}
\value{
A text file containing row-wise activation rules per network component.}
\references{
\emph{DDEPN}\cr
Bender et. al. 2010: Dynamic deterministic effects propagation networks: learning 
 signalling pathways from longitudinal protein array data; Bioinformatics,
 Vol. 26(18), pp. i596-i602\cr
 
\emph{BoolNet}\cr
S. A. Kauffman (1969), Metabolic stability and epigenesis in randomly constructed nets. J. Theor. Biol. 22:437--467.\cr

S. A. Kauffman (1993), The Origins of Order. Oxford University Press.\cr
}
\author{
Silvia von der Heyde
}
\note{
The format of the output file is especially adapted to the \cite{loadNetwork} function format of the \cite{BoolNet-package}.}

\seealso{
\code{\link{create_signetwork}}, \cite{loadNetwork}
}
\examples{
\dontrun{
library(ddepn)
library(BoolNet)

# create example adjacency matrix
example.mat<-matrix(sample(0:2, size=16, replace=TRUE), nrow=4, ncol=4, dimnames=list(x=paste("protein",letters[1:4],sep="_"), y=paste("protein",letters[1:4],sep="_")))

# define output file
example.out<-"exampleOutput.txt"

# convert adjacency matrix to logical rules
adjacencyMatrix_to_logicalRules(adjMatrix=example.mat, outfile=example.out)
}
}

