\name{internal.functions}
\Rdversion{1.1}
\alias{internal.functions}
\alias{get.data}
\alias{plotmatrix}
\alias{get.partitions}
\alias{get.gammashifted}
\alias{transitive.reduction}
\alias{transitive.closure}
\alias{tp.median}
\alias{detailed.to.simple.regulations}
\alias{reverse.direction}

\title{
internal.functions
get.data
plotmatrix
get.gammashifted
get.partitions
transitive.reduction
transitive.closure
tp.median
detailed.to.simple.regulations
reverse.direction
}
\description{
Some internal functions, usually not to be called directly.
}
\usage{
get.data(gammax,mu.bg=0, sd.bg=0.1,
		mu.signal.a=1, sd.signal.a=0.5,
		mu.signal.i=-1, sd.signal.i=0.5,
		stimulus=NULL,TT=10,R.t=4,R.b=3)

plotmatrix(mat,name="")

get.gammashifted(gammax, partitions)

get.partitions(TT, npairs)

transitive.reduction(g)

transitive.closure(g, mat = FALSE, loops = TRUE)

tp.median(dat)

detailed.to.simple.regulations(phi)

reverse.direction(phi, i, switchtype=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
Functions are usually not called directly.
\item{gammax}{todo}
\item{mu.bg, sd.bg, mu.signal.a, sd.signal.a, mu.signal.i, sd.signal.i}{todo}
\item{stimulus}{todo}
\item{R.t,R.b}{todo}
\item{partitions}{todo}
\item{mat}{todo}
\item{name}{todo}
\item{TT}{todo}
\item{npairs}{todo}
\item{g}{todo}
\item{loops}{todo}
\item{dat}{todo}
\item{phi}{todo}
\item{i}{todo}
\item{switchtype}{todo}
 }
\details{
Some internal functions, usually not to be called directly.
The code for \code{transitive.reduction} and \code{transitive.closure} was written 
by Holger Froehlich and is published in the R-package 'nem' available on 
bioconductor. However, to avoid forcing to load the rather big package 'nem' 
the code is included in this package, too.
}
\value{
none
}
\author{
Christian Bender
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

\keyword{ internal }
%\keyword{ ~matrix }% __ONLY ONE__ keyword per line
