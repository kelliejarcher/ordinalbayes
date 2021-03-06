% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.ordinalbayes}
\alias{plot.ordinalbayes}
\title{Trace Plot and/or Density Plot of MCMC Output.}
\usage{
\method{plot}{ordinalbayes}(x, trace = TRUE, density = FALSE, ...)
}
\arguments{
\item{x}{an \code{ordinalBayes} object.}

\item{trace}{a logical value. If TRUE, trace plots are produced for each variable in the chain.}

\item{density}{a logical value. If TRUE, density plots are produced for each variable in the chain.}

\item{...}{other arguments.}
}
\value{
No returned value, called for side effects
}
\description{
Trace Plot and/or Density Plot of MCMC Output.
}
\examples{
\donttest{
data("cesc")
fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],
     model="regressvi",gamma.ind="fixed",
     pi.fixed=0.99, adaptSteps=1000, burnInSteps=1000, nChains=2,
     numSavedSteps=2000, thinSteps=2, seed=26)
 plot(fit)
}
}
