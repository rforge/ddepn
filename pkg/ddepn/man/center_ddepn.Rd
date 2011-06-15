\name{center_ddepn}
\alias{center_ddepn}

\title{
Perform centering of single replicate time courses.
}
\description{
For each experiment in a data matrix (i.e. stimulation/inhibition) subtract the median of 
the curve and shift up to the median of the whole experiment, i.e. median of all biological
and technical replicates.
}
\usage{
center_ddepn(dat)
}
\arguments{
  \item{dat}{The data matrix obtained in \code{\link{hcc1954raw}}, after fcf normalisation. 
    See examples for how to fcf normalise and center the data.}
}
\details{
Used to center the replicate time courses. The time curves should lie as close together as
possible, so that the estimated active/passive Gaussians estimated in \code{\link{ddepn}} are
not obfuscated by the variance between the replicates. FCF normalising and centering will 
yield the data matrix obtained in \code{\link{hcc1954}}. The function is provided for
convenience. It should be traceable, how the matrices were generated.
}
\value{
The centered data matrix.
}
\author{
Christian Bender
}

\seealso{
\code{\link{ddepn}}
\code{\link{hcc1954}}
\code{\link{hcc1954raw}}
}
\examples{
\dontrun{
library(ddepn)
data(hcc1954raw)

## perform FCF normalisation
datfcf <- hcc1954raw / hcc1954fcf * hcc1954fcfmedian

## center the data
datc <- center_ddepn(datfcf)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
