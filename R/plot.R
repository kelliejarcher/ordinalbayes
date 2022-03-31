#' Trace Plot and/or Density Plot of MCMC Output.
#'
#' @param x an \code{ordinalBayes} object.
#' @param trace a logical value. If TRUE, trace plots are produced for each variable in the chain.
#' @param density a logical value. If TRUE, density plots are produced for each variable in the chain.
#' @param ... other arguments.
#'
#' @export
#'
#' @return No returned value, called for side effects
#'
#' @examples
#' \donttest{
#' data("cesc")
#' fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],
#'      model="regressvi",gamma.ind="fixed",
#'      pi.fixed=0.99, adaptSteps=1000, burnInSteps=1000, nChains=2,
#'      numSavedSteps=2000, thinSteps=2, seed=26)
#'  plot(fit)
#' }
plot.ordinalbayes<-function(x, trace = TRUE, density = FALSE, ...) {
  plot(x$results$mcmc, density=density, trace=trace, ...)
}
