\name{create_signetwork}
\alias{create_signetwork}
\alias{perform_network_tests}
\alias{perform_edge_test}
\alias{draw_segments}

\title{
inhibMCMC Analysis: find significant edges from independent MCMC chains
}
\description{
Uses statistical tests (Wilcox Rank-Sum or Ttest) to find significant edges from a number
of independent MCMC runs (performed by passing \code{inference="mcmc"} to the \code{\link{ddepn}}
function call.
}
\usage{
create_signetwork(ret, alpha = 0.001, adj.method = "BY", plot = FALSE, type = "wilcox", alternative = "one.sided", paired = FALSE, ord = NULL)

perform_network_tests(ret, alpha=0.05, plot=FALSE, type="wilcox", alternative="one.sided", paired=FALSE, ord=NULL)

perform_edge_test(ret, from, to, alpha=0.05, plot=TRUE, type="wilcox", alternative="one.sided", paired=FALSE)

draw_segments(tx, rng, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ret}{List. Either the returned object of the \code{\link{ddepn}} call using the inhibMCMC inference
    (argument \code{inference="mcmc"}), or the sublist \code{samplings} from this object. 
  }
  \item{alpha}{
    Numeric in \eqn{[0;1]}. Significance level for the tests.
  }
  \item{adj.method}{
    String. Adjustment for multiple testing. See function \code{p.adjust} from the \emph{stats} package for
    p-value adjustment methods. Most common settings are \code{"none"}, \code{"BH"} or \code{"none"}.  
  }
  \item{plot}{
    Boolean. Should boxplots of the activation/inhibition edge counts (respectively confidences) be plotted?
  }
  \item{type}{
    String. Which test to use, either \code{"wilcox"} or \code{"ttest"}.
  }
  \item{alternative}{
    String. Alternative to test. One out of \code{"one.sided", "two.sided"}.
  }
  \item{paired}{
    Boolean. Paired test? Should normally be unpaired, i.e. FALSE.
  }
  \item{ord}{
    Vector of strings. Defines a node ordering.
  }
  \item{from}{
    String. Defines the source node, from which an edge originates. Used in \code{perform_edge_test} and 
    must be one of the node names of the components in the network.
  }
  \item{to}{
    String. Defines the target node, to which an edge points. Used in \code{perform_edge_test} and 
    must be one of the node names of the components in the network.
   }
  \item{tx}{
    Numeric Vector. Holds the pvalue for alternative greater and alternative less.
  }
  \item{rng}{
    The range of the values that were tested against each other.
  }
  
}

\details{
The main function is \code{create_signetwork}.\cr
Using multiple inhibMCMC chains (a number \emph{L}), a significance testing procedure is used to calculate signfificantly occurring edges.
For each edge, the number of sampled activations and inhibitions is counted in each run, divided by the total number 
of sampled edges, resulting in a 'confidence' for each edge for each run. To find significant edges, a statistical 
test is performed (usually Wilcoxon Rank Sum test) to test the NULL that the means of the \emph{L} activation confidences (c_a)
equal the mean of the \emph{L} inhibition confidences (c_i). The alternatives \eqn{c_a > c_i} and \eqn{c_a < c_i} are 
used to determine the type of interaction.\cr\cr

The functions \code{perform_network_tests}, \code{perform_edge_test} and \code{draw_segments} are usually not called directly.
}
\value{
An adjacency matrix containing the significant edges.
}
\note{
Note that create_signetwork can only be performed for the results of the MCMC inference types and also only if multiple
chains were run. 
}
\references{
Bender 2011, Systematic analysis of time resolved high-throughput data using stochastic network inference methods, 
Dissertation, University of Heidelberg, Faculty of Biology \cr


%% ~put references to the literature/web site here ~
}
\author{
Christian Bender
}


\seealso{
\code{\link{ddepn}}, \code{\link{mcmc_ddepn}}
}
\examples{
\dontrun{
## load package
library(ddepn)
library(multicore)

## sample a network and data
set.seed(1234)
n <- 6
signet <- signalnetwork(n=n, nstim=2, cstim=0, prop.inh=0.2)
phit <- signet$phi
stimuli <- signet$stimuli
dataset <- makedata(phit, stimuli, mu.bg=1200, sd.bg=400, mu.signal.a=2000, sd.signal.a=1000)

## use original network as prior matrix
## reset all entries for inhibiting edges 
## to -1
B <- phit
B[B==2] <- -1

## perform inhibMCMC inference, using 4 CPU cores to get 4 MCMC chains 
ret <- ddepn(dataset$datx, phiorig=phit, maxiterations=300, burnin=50,
            plotresults=FALSE, inference="mcmc",
            usebics=FALSE, priortype="laplaceinhib", lambda=0.01, B=B,
            multicores=TRUE, cores=4)
            
## get the signal network containing the significant edges
net <- create_signetwork(ret, alpha=0.05, adj.method="none", plot=FALSE, type="wilcox", alternative="one.sided", paired=FALSE, ord=NULL)

## compare inferred net and prior
par(mfrow=c(1,2))
plotdetailed(net, main="inferred", stimuli=ret$samplings[[1]]$stimuli )
plotdetailed(ret$samplings[[1]]$phi.orig, main="prior", stimuli=ret$samplings[[1]]$stimuli )



}
}
