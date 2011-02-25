\name{format_ddepn}
\alias{format_ddepn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs renaming of the colnames of dataset hcc1954 to make it 
suitable as input for the inference.
}
\description{
A little helper that strips of the information on biological replicates in the 
example data matrix \code{\link{hcc1954}}. After stripping the matrix can be directly used
as input for \code{\link{ddepn}}.
}
\usage{
format_ddepn(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data matrix in the format described in \code{\link{ddepn}}, except for
    the column names, which still contain a label for the biological replicate. See 
    the details for an example.
    
%%     ~~Describe \code{dat} here~~
}
}
\details{
  Example for column names holding information on biological replicates:
  \tabular{lcccccccc}{
    \tab EGF-3_1\tab EGF-5_1 \tab EGF-3_1 \tab EGF-5_1 \tab EGF-3_2 \tab EGF-5_2 \tab EGF-3_2 \tab EGF-5_2 \cr
    EGF ...\tab \tab  \tab  \tab  \tab  \tab  \tab  \tab \cr
    X   ...\tab \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
    AKT ...\tab \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
  }
  The format is \code{STIMULUS-biolreplicate_time}, i.e. \code{EGF-3_1} means experiment with
  EGF stimulation, biological replicate labeled as 3 at time point 1. The biological replicate
  should be labeled by numbers only.
}
\author{
Christian Bender
}

\seealso{
\code{\link{hcc1954}}
}
\examples{
  \dontrun{
    ## load package
    library(ddepn)
    data(hcc1954)
    format_ddepn(hcc1954)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
