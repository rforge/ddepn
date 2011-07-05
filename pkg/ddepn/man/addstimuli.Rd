\name{addstimuli}
\alias{addstimuli}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add dummy stimulus nodes to the data matrix.
}
\description{
Add the treatments (stimuli, inhibitions) defined in the column names of the data matrix to the data matrix as dummy nodes, if not there.
}
\usage{
addstimuli(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  The data matrix. See \code{\link{ddepn}} for more information.
}
}
\details{
Uses the information from the column names of the data matrix to derive the input treatments 
(either stimuli or inhibition) for the experiment. Each stimulus is added as node to the 
network to be inferred, and thus as row to the data matrix. For example, activating treatments
could be stimulation of cells by growth factors, inhibiting treatments could be treatment 
with a specific drug. 
}
\value{
The modified data matrix.
}
\author{
Christian Bender
}

\seealso{
\code{\link{ddepn}}
}
\examples{
\dontrun{
TODO
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ internal }
