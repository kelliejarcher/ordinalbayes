#' Extract Model Coefficients
#'
#' @param object an \code{ordinalbayes} object.
#' @param method The default is \code{method=mean} which estimates the mean of each parameter in the MCMC chain. Other options are \code{method=median} or any other relevant summary function.
#' @param ... other arguements.
#'
#' @return \item{alpha }{Summary estimates for the thresholds}
#' @return \item{zeta }{Summary estimates for the unpenalized covariates. Only available if unpenalized covariates were included in the fitted model.}
#' @return \item{beta}{Summary estimates for the penalized covariates}
#' @return \item{gamma}{Summary estimates for the variable inclusion indicators. Not available when \code{model="lasso"}}
#'
#' @seealso \code{\link{ordinalbayes}}, \code{\link{print.ordinalbayes}}, \code{\link{summary.ordinalbayes}}, \code{\link{predict.ordinalbayes}}
#'
#' @examples
#' \donttest{
#' data("cesc")
#' fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],
#'          model="regressvi", gamma.ind="fixed", pi.fixed=0.99,
#'          adaptSteps=1000, burnInSteps=1000, nChains=2,
#'          numSavedSteps=2000, thinSteps=2, seed=26)
#' coef(fit)
#' }
#' @keywords methods
#' @export
#' @method coef ordinalbayes
coef.ordinalbayes <-
  function(object, method=mean, ...) {
    summary <- summary(object)
    alpha <-apply(summary$alphamatrix,2,function(x) eval(expression(method(x))))
    names(alpha)<-paste("alpha",1:length(alpha),sep="")
    if (dim(object$w)[2]!=0) {
      zeta <-apply(summary$zetamatrix,2,function(x) eval(expression(method(x))))
      names(zeta)<-dimnames(object$w)[[2]]
    }
    if (object$model=="lasso") {
      beta<-apply(summary$betamatrix,2,function(x) eval(expression(method(x))))
      names(beta)<-object$featureNames
      if (dim(object$w)[2]==0) {
        output<-list(alpha, beta)
        names(output)<-c("alpha", "beta")
      } else {
        output<-list(alpha, zeta, beta)
        names(output)<-c("alpha", "zeta", "beta")
      }
    } else {
      beta<-apply(summary$betamatrix,2,function(x) eval(expression(method(x))))
      names(beta)<-object$featureNames
      gamma<-apply(summary$gammamatrix,2,function(x) eval(expression(method(x))))
      names(gamma)<-object$featureNames
      if (dim(object$w)[2]==0) {
        output<-list(alpha, beta, gamma)
        names(output)<-c("alpha", "beta", "gamma")
      } else {
        output<-list(alpha, zeta, beta, gamma)
        names(output)<-c("alpha", "zeta", "beta", "gamma")
      }
    }
    output
  }
