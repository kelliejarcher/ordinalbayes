#' Summarize an Ordinal Bayes Object.
#'
#' @param object A fitted \code{ordinalbayes} object.
#' @param epsilon   a small positive value that is close to 0 for testing an interval null hypothesis for the beta parameters.
#' @param ... other arguments.
#'
#' @return \item{alphamatrix}{The MCMC output for the threshold parameters.}
#' @return \item{betamatrix}{The MCMC output for the penalized parameters.}
#' @return \item{zetamatrix}{The MCMC output for the unpenalized parameters (if included).}
#' @return \item{gammamatrix}{The MCMC output for the variable inclusion parameters (not available for lasso).}
#' @return \item{gammamean}{The posterior mean of the variable inclusion indicators (not available for lasso)for the variable inclusion indicators (not available for lasso).}
#' @return \item{gamma.BayesFactor}{Bayes factor for the variable inclusion indicators (not available for lasso).}
#' @return \item{Beta.BayesFactor}{Bayes factor for the penalized parameters where the interval null is tested using epsilon.}
#' @return \item{lambdamatrix}{The MCMC output for the penalty parameter (not available for normalss).}
#'
#' @export
#'
#' @seealso \code{\link{ordinalbayes}}, \code{\link{print.ordinalbayes}}, \code{\link{coef.ordinalbayes}}, \code{\link{predict.ordinalbayes}}
#'
#' @examples
#' \donttest{
#' data("cesc")
#' fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],
#'          model="regressvi", gamma.ind="fixed", pi.fixed=0.99,
#'          adaptSteps=1000, burnInSteps=1000, nChains=2,
#'          numSavedSteps=2000, thinSteps=2, seed=26)
#' summary.fit<-summary(fit)
#' names(summary.fit)
#' names(which(summary.fit$Beta.BayesFactor>5))
#' names(which(summary.fit$gamma.BayesFactor>5))
#' }
#' @importFrom stats pnorm
#' @method summary ordinalbayes
summary.ordinalbayes <-
  function(object, epsilon=0.1, ...) {
    a<-object$a
    b<-object$b
    model<-object$model
    if (length(object$results$mcmc)>1) {
      mcmcChain <-object$results$mcmc
      for(i in 2:length(object$results$mcmc)) {
        mcmcChain[[i]] <-rbind(mcmcChain[[i-1]], object$results$mcmc[[i]])
      }
      mcmcChain <- mcmcChain[[length(object$results$mcmc)]]
    } else {
      mcmcChain <-as.matrix(object$results$mcmc)
    }
    if (dim(object$w)[2]!=0) {
      w<-as.matrix(object$w,drop=FALSE)
      varNames<-dimnames(w)[[2]]
      zetamatrix <- mcmcChain[,match(varNames, dimnames(mcmcChain)[[2]]),drop=FALSE]
      colnames(zetamatrix)<-varNames
    }
    featureNames<-object$featureNames
    if (model!="normalss") {
      lambdamatrix <- mcmcChain[,grep("lambda", dimnames(mcmcChain)[[2]])]
    }
    alphamatrix <- mcmcChain[,grep("alpha", dimnames(mcmcChain)[[2]])]
    if (model=="lasso") {
      betamatrix <- mcmcChain[,grep("beta",dimnames(mcmcChain)[[2]])]
      PriorOdds2 = b^a/((b+epsilon)^a-b^a)
    } else if (!is.na(match(model,c("normalss","dess","regressvi")))) {
      gamma.ind<-object$gamma.ind
      if (c("fixed","random")[charmatch(gamma.ind,c("fixed","random"))]=="fixed") {
        p <- object$pi.fixed
      } else {
        p = beta(object$c.gamma+1,object$d.gamma)/beta(object$c.gamma,object$d.gamma)
      }
      if (model=="regressvi") {
        betamatrix <- mcmcChain[,grep("beta",dimnames(mcmcChain)[[2]])]
        PriorOdds2 = (p*b^a*gamma(a+1))/(p*gamma(a+1)*(b+epsilon)^a-p*gamma(a+1)*b^a+(1-p)*a*gamma(a)*(b+epsilon)^a)
      }
      if (model=="normalss") {
        sigma0 = sqrt(object$sigma2.0) # sigma0 and sigma1 are standard deviations for the two normal distributions constitute the priors for beta
        sigma1 = sqrt(object$sigma2.1)
        PriorOdds2 <- (2*(1-p)*(1-pnorm(epsilon/sigma0))+2*p*(1-pnorm(epsilon/sigma1)))/((1-p)*(pnorm(epsilon/sigma0)-pnorm(-epsilon/sigma0))+p*(pnorm(epsilon/sigma1)-pnorm(-epsilon/sigma1)))
        betamatrix <- mcmcChain[,grep("beta",dimnames(mcmcChain)[[2]])]
      }
      if (model=="dess") {
        PriorOdds2 <- (p*b^a*gamma(a+1)+(1-p)*a*gamma(a)*(b+epsilon)^a*exp(-object$lambda0*epsilon))/(p*gamma(a+1)*((b+epsilon)^a-b^a)+(1-p)*a*gamma(a)*(1-exp(-object$lambda0*epsilon))*(b+epsilon)^a)
        betamatrix <- mcmcChain[,grep("beta",dimnames(mcmcChain)[[2]])]
      }
      gammamatrix <- mcmcChain[,grep("gamma",dimnames(mcmcChain)[[2]])]
      ### select variable using top N gamma means
      gammamean = apply(gammamatrix, 2,mean)
      names(gammamean) <- object$featureNames
      ### estimate Bayes factor for gamma
      gammasum = apply(gammamatrix, 2,sum)
      PostOdds = gammasum/(dim(gammamatrix)[1]-gammasum)
      PriorOdds = p/(1-p)
      gamma.BayesFactor<-PostOdds/PriorOdds
    }
    ### estimate Bayes factor for beta
    BetgmaPost1 = apply(betamatrix, 2, function(x) sum(abs(x) > epsilon))
    BetgmaPost2 = apply(betamatrix, 2, function(x) sum(abs(x) <= epsilon))
    PostOdds2 = BetgmaPost1/BetgmaPost2
    Beta.BayesFactor = PostOdds2/PriorOdds2
    names(Beta.BayesFactor)<-object$featureNames
    if (model=="lasso") {
      if (dim(object$w)[2]==0) {
        output<-list(alphamatrix, betamatrix, Beta.BayesFactor, lambdamatrix)
        names(output)<-c("alphamatrix", "betamatrix", "Beta.BayesFactor", "lambdamatrix")
      } else {
        output<-list(alphamatrix, betamatrix, zetamatrix, Beta.BayesFactor, lambdamatrix)
        names(output)<-c("alphamatrix", "betamatrix", "zetamatrix", "Beta.BayesFactor", "lambdamatrix")
      }
    } else {
      names(gamma.BayesFactor)<-object$featureNames
      if (dim(object$w)[2]==0) {
        if (model=="normalss") {
          output<-list(alphamatrix, betamatrix, gammamatrix, gammamean, gamma.BayesFactor, Beta.BayesFactor)
          names(output)<-c("alphamatrix", "betamatrix", "gammamatrix", "gammamean","gamma.BayesFactor","Beta.BayesFactor")
        } else {
          output<-list(alphamatrix, betamatrix, gammamatrix, gammamean, gamma.BayesFactor, Beta.BayesFactor, lambdamatrix)
          names(output)<-c("alphamatrix", "betamatrix", "gammamatrix", "gammamean","gamma.BayesFactor","Beta.BayesFactor", "lambdamatrix")
        }
      } else {
        if (model=="normalss") {
          output<-list(alphamatrix, betamatrix, zetamatrix, gammamatrix, gammamean, gamma.BayesFactor, Beta.BayesFactor)
          names(output)<-c("alphamatrix", "betamatrix", "zetamatrix", "gammamatrix", "gammamean","gamma.BayesFactor","Beta.BayesFactor")
        } else {
          output<-list(alphamatrix, betamatrix, zetamatrix, gammamatrix, gammamean, gamma.BayesFactor, Beta.BayesFactor, lambdamatrix)
          names(output)<-c("alphamatrix", "betamatrix", "zetamatrix", "gammamatrix", "gammamean","gamma.BayesFactor","Beta.BayesFactor", "lambdamatrix")
        }
      }
    }
    output
  }
