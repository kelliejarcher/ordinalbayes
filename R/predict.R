#' Predicted Probabilities and Class for an Ordinal Bayes Fit.
#'
#' @aliases fitted.ordinalbayes
#'
#' @param object an \code{ordinalBayes} fitted object.
#' @param neww an optional formula that includes the unpenalized variables to use for predicting the response. If omitted, the training data are used.
#' @param newdata an optional data.frame that minimally includes the unpenalized variables to use for predicting the response. If omitted, the training data are used.
#' @param newx an optional matrix of penalized variables to use for predicting the response. If omitted, the training data are used.
#' @param model.select when \code{"average"} (default) is used, the mean coefficient values over the MCMC chain are used to estimate fitted probabilities; when  \code{"median"} is used, the median coefficient values over the MCMC chain are used to estimate fitted probabilities;  when \code{"max.predicted.class"} is used, each step in the chain is used to calculate fitted probabilities and the class. The predicted class is that attaining the maximum fitted probability.
#' @param ... other arguments.
#'
#' @return \item{predicted}{a matrix of predicted probabilities from the fitted model.}
#' @return \item{class}{a vector containing the predicted class taken as that class having the largest predicted probability.}

#' @export
#'
#' @seealso \code{\link{ordinalbayes}}, \code{\link{coef.ordinalbayes}}, \code{\link{summary.ordinalbayes}}, \code{\link{print.ordinalbayes}}
#'
#' @examples
#' \donttest{
#' data("cesc")
#' fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],,
#'      model="regressvi", gamma.ind="fixed", pi.fixed=0.99,
#'      adaptSteps=1000, burnInSteps=1000, nChains=2, numSavedSteps=2000,
#'      thinSteps=2, seed=26)
#' phat<-predict(fit)
#' table(phat$class, cesc$Stage)
#' }
#' @importFrom stats median model.frame
predict.ordinalbayes <-
  function(object, neww = NULL, newdata, newx = NULL, model.select = "average", ...) {
    Y <- object$y
    w <- object$w
    x <- object$x
    levels <- levels(Y)
    if (length(object$results$mcmc)>1) {
      mcmcChain <-object$results$mcmc
      for(i in 2:length(object$results$mcmc)) {
        mcmcChain[[i]] <-rbind(mcmcChain[[i-1]], object$results$mcmc[[i]])
      }
      mcmcChain <- mcmcChain[[length(object$results$mcmc)]]
    } else {
      mcmcChain <-as.matrix(object$results$mcmc)
    }
    featureNames<-object$featureNames
    betamatrix <- mcmcChain[,grep("bet",dimnames(mcmcChain)[[2]])]
    alphamatrix <- mcmcChain[,grep("alpha",dimnames(mcmcChain)[[2]])]
    if (dim(w)[2] != 0) {
      w<-as.matrix(object$w,drop=FALSE)
      varNames<-dimnames(w)[[2]]
      zetamatrix <- mcmcChain[,match(varNames, dimnames(mcmcChain)[[2]]),drop=FALSE]
      colnames(zetamatrix)<-varNames
    }
    k <- length(unique(Y))
    if (!is.null(neww))
      if (neww == ~1) {
        m <- model.frame(neww)
      } else {
        m <- model.frame(neww, newdata)
        neww <- model.matrix(neww, m)
        neww <- neww[, -grep("(Intercept)", dimnames(neww)[[2]]),
                     drop = FALSE]
      }
    if (!is.null(newx))
      newx <- as.matrix(newx)
    if (is.null(neww) & is.null(newx)) {
      neww <- object$w
      newx <- object$x
    }
    n <- max(dim(neww)[1], dim(newx)[1])
    if (!is.null(newx) & identical(newx, x)) {
      if (object$scale) {
        sd <- apply(newx, 2, sd)
        for (i in 1:dim(newx)[2]) {
          if (sd[i] == 0) {
            newx[, i] <- scale(newx[, i], center = TRUE,
                               scale = FALSE)
          }
          else {
            newx[, i] <- scale(newx[, i], center = TRUE,
                               scale = TRUE)
          }
        }
      }
    } else if (!is.null(newx) && object$scale) {
      newx <- rbind(x, newx)
      sd <- apply(newx, 2, sd)
      for (i in 1:dim(newx)[2]) {
        if (sd[i] == 0) {
          newx[, i] <- scale(newx[, i], center = TRUE,
                             scale = FALSE)
        }
        else {
          newx[, i] <- scale(newx[, i], center = TRUE,
                             scale = TRUE)
        }
      }
      newx <- matrix(newx[-(1:dim(x)[1]), ], ncol = dim(x)[2])
    }
    if (model.select=="average") {
      coef<-coef(object, method=mean)
      alpha<-coef$alpha
      beta<-coef$beta
      if (dim(w)[2]!=0)
        zeta<-coef$zeta
    } else if (model.select=="median") {
      coef<-coef(object, method=median)
      alpha<-coef$alpha
      beta<-coef$beta
      if (dim(w)[2]!=0)
        zeta<-coef$zeta
    }
    if (model.select!="max.predicted.class") {
      if (dim(w)[2] != 0) {
        if (is.null(x)) {
          Xb <- neww %*% zeta
        } else if (!is.null(x)) {
          Xb <- neww %*% zeta + newx %*% beta
        }
      } else if (!is.null(x)) {
        Xb <- newx %*% beta
      } else {
        Xb <- 0
      }
      z <- matrix(ncol = k - 1, nrow = n)
      for (i in 1:(k - 1)) {
        z[, i] <- alpha[i] - Xb
      }
      pi <- matrix(ncol = k, nrow = n)
      for (i in 1:k) {
        if (i == 1) {
          pi[, i] <- exp(z[, i])/(1 + exp(z[, i]))
        }  else if (i <= k - 1) {
          pi[, i] <- exp(z[, i])/(1 + exp(z[, i])) -
            exp(z[, i - 1])/(1 + exp(z[, i - 1]))
        }  else if (i == k) {
          pi[, i] <- 1 - exp(z[, i - 1])/(1 + exp(z[,
                                                    i - 1]))
        }
      }
      class <- levels[apply(pi, 1, which.max)]
      output <- list(predicted = pi, class = class)
    } else {
      if (dim(w)[2] != 0) {
        if (is.null(x)) {
          XBeta <- zetamatrix%*%t(neww)
        } else if (!is.null(x)) {
          XBeta <- zetamatrix%*%t(neww)+betamatrix%*%t(newx)
        }
      } else if (!is.null(x)) {
        XBeta<-betamatrix%*%t(newx)
      } else {
        XBeta<-0
      }
      z <- array(dim=c(k - 1, dim(betamatrix)[1], n))
      for (i in 1:(k-1)) {
        z[i,,] <- alphamatrix[,i] - XBeta
      }
      pi <- array(dim=c(k, dim(betamatrix)[1], n))
      for (i in 1:k) {
        if (i == 1) {
          pi[i,,] <- exp(z[i,,])/(1 + exp(z[i,,]))
        }  else if (i <= k - 1) {
          pi[i,,] <- exp(z[i,,])/(1 + exp(z[i,,])) -
            exp(z[i - 1,,])/(1 + exp(z[i - 1,,]))
        }  else if (i == k) {
          pi[i,,] <- 1 - exp(z[i - 1,,])/(1 + exp(z[
            i - 1,,]))
        }
      }
      class.temp <- matrix(nrow=n, ncol=dim(betamatrix)[1])
      class<-numeric()
      for (j in 1:n){
        class.temp[j,] <- apply(pi[,,j],2,which.max)
        table.temp<-table(class.temp[j,])
        class.no<-as.numeric(names(table.temp)[which.max(table.temp)])
        class[j]<-levels[class.no]
      }
      output <- list(class=class)
    }
    output
  }
