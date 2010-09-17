\name{makedata}
\alias{makedata}
\title{
Generate artificial dataset given a network and list of stimuli.
}
\description{
Dataset generation using an input network, list of stimuli and various experiment parameters.
}
\usage{
makedata(phi, stimuli, R.t = 3, R.b = 3, TT = 10, mu.bg = 0, sd.bg = 0.1, mu.signal.a = 2, sd.signal.a = 0.5, mu.signal.i = -2, sd.signal.i = 0.5)
}
\arguments{
  \item{phi}{Adjacency matrix.The input network.
  }
  \item{stimuli}{List containing all stimuli.
  }
  \item{R.t}{Integer. Number of technical replicates. So far, no distinction is done between technical and 
    biological replicates, so in total a number of \eqn{R.t \cdot R.b} replicates is generated. 
  }
  \item{R.b}{Integer. Number of biological replicates. So far, no distinction is done between technical and 
    biological replicates, so in total a number of \eqn{R.t \cdot R.b} replicates is generated.
  }
  \item{TT}{Integer.Number of timepoints in the experiment.
  }
  \item{mu.bg}{Double. Mean background intensity (for passive state of the protein.
  }
  \item{sd.bg}{Double. Sd for background intensity.
  }
  \item{mu.signal.a}{Double. Mean intensity for activation.
  }
  \item{sd.signal.a}{Double. Sd for activation intensity.
  }
  \item{mu.signal.i}{Double. Mean intensity for inhibition.
  }
  \item{sd.signal.i}{Double. Sd for inhibition intensity.
  }
}
\details{
Generates a dataset from the given network, stimuli and parameters.
}
\value{
A list containing the dataset:
	\item{datx}{The data matrix.}
	\item{gammax}{Matrix of true state transitions.}
	\item{stimuli}{The list of stimuli.}
	\item{phi}{The network.}
	\item{R.t}{Technical replicates.}
	\item{R.b}{Biological replicates.}
	\item{TT}{Number of timepoints.} 
}
\author{
Christian Bender
}
\examples{
\dontrun{
library(ddepn)
dataset <- makedata(matrix(sample(c(0,1,2),9,replace=TRUE),nrow=3,dimnames=list(LETTERS[1:3],LETTERS[1:3])), list(A=1))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ internal }
