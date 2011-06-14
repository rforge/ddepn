\name{bestgam}
\alias{bestgam}
\alias{esub}
\title{
Fit an optimal smoothing spline to a data vector y, indexed by time vector tp.
}
\description{
For a data vector \code{y}, where \code{tp} assigns a time point to each value in \code{y}, a smoothing spline
is fitted through the data points. To find an optimal spline, a series of fits is performed with increasing number 
of degrees of freedom and the optimal fit is selected via ANOVA testing and selection using the \code{selection.criterion}.
}
\usage{
bestgam(y, tp, xn, selection.criterion = "aic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    Numeric vector. The data points.
}
  \item{tp}{
    Numeric vector. The corresponding time points.
}
  \item{xn}{
    Numeric vector. The time values for the prediction.
}
  \item{selection.criterion}{
    String. Either \code{"aic"} for Akaike's Information Criterion, \code{"BIC"} for Bayesian Information Criterion 
    or \code{"pvalue"} for model selection according to the p-value.
}
}
\details{
Several smoothing splines with increasing number of degrees of freedom are fit to the data (maximal number of df is 
one less than time points). Additionally, a constant and linear model fit is produced. The linear and smoothing spline 
fits are each compared to the constant fit via analysis of variance/deviance (ANOVA) to decide whether significant
improvement of the complex model fits compared to the simple constant fit were got. This can also be seen as 
identification of a significant effect over time. To select a specific curve, the resulting p-values from the ANOVA are
used. Either directly by selecting the minimal p-value of all fit comparisons, or by plugging in the p-value into AIC or 
BIC and maximising the respective statistic.\cr

After the optimal fit is selected, a prediciton of the time profile is calculated, using the time vector \code{xn} as 
predictor variable in the fit. \cr

\code{esub} is a helper for the gam fitting, usually not called directly.
}
\value{
Returns a vector with the predicted values for the time vector \code{xn}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Christian Bender\cr
Tim Beissbarth\cr
Stephan Gade\cr
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
## load package
library(ddepn)
library(multicore)
library(gam)

## sample a network and data
set.seed(1234)
n <- 6
signet <- signalnetwork(n=n, nstim=2, cstim=0, prop.inh=0.2)
phit <- signet$phi
stimuli <- signet$stimuli
dataset <- makedata(phit, stimuli, mu.bg=1200, sd.bg=400, mu.signal.a=2000, sd.signal.a=1000)

## select an experiment

prepare_bestgam <- function(dataset, experiment, protein) {
  stimuli <- dataset$stimuli
  dat <- dataset$datx
  stim <- names(stimuli[[experiment]])
  sel <- grep(paste("^",stim,"_.*$",sep=""), colnames(dat))
  y <- dat[protein,sel]
  tp <- as.numeric(sub("^.*_","",colnames(dat)[sel]))
  #ord <- order(tp)
  #y <- y[ord]
  #tp <- tp[ord]
  xn <- seq(0,max(tp),0.1)
  return(list(y=y, tp=tp, xn=xn))
}

## no effect observed in experiment 1, protein 2
lst <- prepare_bestgam(dataset, 1, 2)
yhat <- bestgam(lst$y, lst$tp, lst$xn, "aic")
plot(lst$y~lst$tp)
lines(lst$xn,yhat)
 
## effect observed in experiment 1, protein 5
lst <- prepare_bestgam(dataset, 1, 5)
yhat <- bestgam(lst$y, lst$tp, lst$xn, "aic")
plot(lst$y~lst$tp)
lines(lst$xn,yhat)
 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
