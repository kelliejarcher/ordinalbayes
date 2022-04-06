#' Print MCMC Summary Statistics
#'
#' @param x A fitted ordinalbayes object.
#' @param ... other arguments.
#'
#' @return Matrix with the summaries from the \code{x$results} object which of class 'runjags'.
#' Columns include Lower95, Median, Upper95, Mean, SD, Mode, MCerr, MC%ofSD, SSeff, AC.20, and psrf
#' @export
#'
#' @seealso \code{\link{ordinalbayes}}, \code{\link{summary.ordinalbayes}}, \code{\link{coef.ordinalbayes}}, \code{\link{predict.ordinalbayes}}
#'
#' @examples
#' \donttest{
#' data("cesc")
#' fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],
#' model="regressvi", gamma.ind="fixed", pi.fixed=0.99, adaptSteps=1000,
#' burnInSteps=1000, nChains=2, numSavedSteps=2000, thinSteps=2)
#' print(fit)
#' }
#' @method print ordinalbayes
print.ordinalbayes<-function(x, ...) {
  res <- x$results$summaries
  print(res)
}
