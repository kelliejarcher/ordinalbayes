#' Ordinal Bayesian Regression Models for High-Dimensional Data
#'
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The left side of the formula is the ordinal outcome while the variables on the right side of the formula are the covariates that are not included in the penalization process. Note that if all variables in the model are to be penalized, an intercept only model formula should be specified.
#' @param data an optional data.frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
#' @param x an optional matrix of predictors that are to be penalized in the model fitting process.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param center logical, if TRUE the penalized predictors are centered.
#' @param scale logical, if TRUE the penalized predictors are scaled.
#' @param a hyperprior for the penalty parameter lambda which is Gamma with parameters \code{a} and \code{b}.
#' @param b hyperprior for the penalty parameter lambda which is Gamma with parameters \code{a} and \code{b}.
#' @param model Specify which penalized ordinal model to fit  as "regressvi", "lasso", "dess", or "normalss".
#' @param gamma.ind indicates whether prior for the variable inclusion indicators is "fixed" or "random" (for models "regressvi", "dess", or "normalss").
#' @param pi.fixed constant prior for the variable inclusion indicators is when \code{gamma.ind="fixed"}.
#' @param c.gamma hyperprior for the variable inclusion indicators is when \code{gamma.ind="random"}.
#' @param d.gamma hyperprior for the variable inclusion indicators is when \code{gamma.ind="random"}.
#' @param alpha.var variance for alpha_k thresholds in the MCMC chain (default 10).
#' @param sigma2.0 variance for the spike when \code{model="normalss"} (set to some small positive value).
#' @param sigma2.1   variance for the slab when \code{model="normalss"} (set to some large positive value).
#' @param coerce.var variance associated with any unpenalized predictors in the MCMC chain (default 10).
#' @param lambda0 parameter value for the spike when \code{model="dess"}.
#' @param nChains number of parallel chains to run (default 3)
#' @param adaptSteps number of iterations for adaptation (default 5,000).
#' @param burnInSteps number of iterations of the Markov chain to run (default 5,000).
#' @param numSavedSteps number of saved steps per chain (default 9,999).
#' @param thinSteps thinning interval for monitors (default 3).
#' @param parallel logical, run the MCMC on multiple processors (default TRUE).
#' @param seed integer, seed to ensure reproducibility.
#' @param quiet logical, when TRUE, suppress output of JAGS (or rjags) when updating models
#'
#' @return \code{results} An object of class runjags
#' @return \code{call} Model call
#' @return \code{model} Name of the ordinal model that was fit
#' @return \code{a} Value the user specified for \code{a}
#' @return \code{b} Value the user specified for \code{b}
#' @return \code{featureNames} Names of the penalized predictors
#' @return \code{center} Value the user specified for \code{center}
#' @return \code{scale} Value the user specified for \code{scale}
#' @return \code{y} Observed ordinal response
#' @return \code{x} Matrix of penalized predictors used in model fitting
#' @return \code{w} Matrix of unpenalized predictors used in model fitting
#' @return \code{gamma.ind} Value the user specified for \code{gamma.ind}
#' @return \code{pi.fixed} Value the user specified for \code{pi.fixed} if \code{gamma.ind="fixed"}
#' @return \code{c.gamma} Value the user specified for \code{c.gamma} if \code{gamma.ind="random"}
#' @return \code{d.gamma} Value the user specified for \code{d.gamma} if \code{gamma.ind="random"}
#' @return \code{sigma2.0} Value the user specified for \code{sigma2.0} if \code{model="normalss"}
#' @return \code{sigma2.1} Value the user specified for \code{sigma2.1} if \code{model="normalss"}
#' @return \code{lambda0} value the user specified for \code{lambda0} if \code{model="dess"}
#' @export
#'
#' @seealso \code{\link{print.ordinalbayes}}, \code{\link{summary.ordinalbayes}}, \code{\link{coef.ordinalbayes}}
#'
#' @examples
#' \donttest{
#' # The number of adaptSteps, burnInSteps, and numSavedSteps was reduced for package testing
#' data("cesc")
#' data(reducedSet)
#' fit<-ordinalbayes(Stage~1, data=cesc, x=cesc[,5:45],
#'          model="regressvi", gamma.ind="fixed", pi.fixed=0.99,
#'          adaptSteps=1000, burnInSteps=1000, nChains=2,
#'          numSavedSteps=2000, thinSteps=2, seed=26)
#' }
#' @keywords models
#' @keywords regression
#' @importFrom stats model.matrix model.response rnorm
ordinalbayes <-
  function(formula, data, x = NULL, subset, center=TRUE, scale=TRUE, a=0.1, b=0.1, model="regressvi", gamma.ind="fixed", pi.fixed=0.05, c.gamma=NULL, d.gamma=NULL, alpha.var=10,
           sigma2.0=NULL, sigma2.1=NULL, coerce.var=10, lambda0=NULL, nChains=3, adaptSteps=5000, burnInSteps=5000,  numSavedSteps=9999, thinSteps=3, parallel=TRUE, seed=NULL, quiet=FALSE) {
    mf <- match.call(expand.dots = FALSE)
    cl <- match.call()
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.response(mf)
    if (is.na(match(model,c("lasso", "normalss", "dess", "regressvi")))) {
      stop("model must be either 'lasso', 'normalss', 'dess' or 'regressvi'\n")
    }
    if (!is.na(match(model,c("normalss","dess","regressvi"))) & is.na(match(c("fixed","random")[charmatch(gamma.ind,c("fixed","random"))],c("fixed","random")))) {
      stop("gamma.ind must be either 'fixed' or 'random' for models 'normalss', 'dess' and 'regressvi'\n")
    }
    if (c("fixed","random")[charmatch(gamma.ind,c("fixed","random"))]=="fixed" & (is.na(pi.fixed) | pi.fixed<=0 | pi.fixed>=1) ) {
      stop("pi.gamma must be a value in (0, 1) when gamma='fixed'\n")
    }
    if (!is.na(match(model,c("normalss","dess","regressvi"))) & c("fixed","random")[charmatch(gamma.ind,c("fixed","random"))]=="random" & (is.null(c.gamma) | is.null(d.gamma))) {
      stop("c.gamma and d.gamma must be specified as parameters of the Beta distribution when gamma.ind is 'random' for the normal spike-and-slab (normalss), double exponential spike-and-slab (dess), and regression-based variable inclusion indicator (regressvi) models\n")
    }
    if (!is.na(match(model,c("normalss"))) && (is.null(sigma2.0) | is.null(sigma2.1))) {
      stop("sigma2.0 and sigma2.1 must be specified for normal spike-and-slab model (normalss)")
    }
    if (!is.na(match(model,c("dess"))) && (is.null(lambda0))) {
      stop("lambda0 must be specified for double exponential spike-and-slab model (dess)")
    }
    if (class(Y)[1]=="ordered") {
      y<-matrix(as.numeric(Y),ncol=1)
      levels<-levels(Y)
    } else {
      stop("Response should be an ordered factor.\n")
    }
    if (quiet) {
      runjags::runjags.options(silent.runjags = TRUE, silent.jags = TRUE)
    } else {runjags::runjags.options(silent.runjags = TRUE, silent.jags = FALSE)}
    nChains <- round(nChains)
    w <- model.matrix(mt, mf)
    if (!is.null(x)) {
      if (missing(subset))
        r <- TRUE
      else {
        e <- substitute(subset)
        r <- eval(e, data)
        if (!is.logical(r))
          stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
      }
      if (class(x)[1] == "character") {
        nl <- as.list(1:ncol(data))
        names(nl) <- names(data)
        vars <- eval(substitute(x), nl, parent.frame())
        x <- data[r, vars, drop = FALSE]
        x <- as.matrix(x)
      }
      else if (class(x)[1] == "matrix" || class(x)[1] == "data.frame") {
        x <- x[r, , drop = FALSE]
        x <- as.matrix(apply(x, 2, as.numeric))
      }
      P<-dim(x)[2]
      if (scale) {
        if (center) {
          x <- scale(x, center = TRUE, scale = TRUE)
        } else {
          x <- scale(x, center = FALSE, scale = TRUE)
        }
      }
    }
    if (dim(y)[1]!=dim(x)[1]) stop("Number of rows in x does not equal the length of the ordinal response vector")
    is.intercept <- grep("Intercept", dimnames(w)[[2]])
    if (length(is.intercept) == 1) {
      w <- w[, -is.intercept, drop = FALSE]
    }
    zeta <- rep(0, dim(w)[2])
    N<-dim(w)[1]
    if (N==0) stop("sample size is 0")
    Q<-dim(w)[2]
    k<-length(levels)
    alpha.precision <- 1/alpha.var
    coerce.precision <- 1/coerce.var
    alpha <- numeric()
    pi.0 <- table(Y)/length(Y)
    tab <- table(Y)
    alpha <- log(cumsum(pi.0)/(1 - cumsum(pi.0)))[1:(k - 1)]
    #nIter <- ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain
    nIter <- numSavedSteps
    inits <- list()
    if (is.null(seed)) {seed<-1}
    set.seed(seed)
    if (model=="lasso") {
      if (dim(w)[2]==0) {
        modelname <- lasso
        file <- runjags::read.jagsfile(lasso)
        jagsdf <- list(y=y, X=x, N=N, P=P, k=k, a=a, b=b, alpha.precision=alpha.precision)
      } else {
        modelname <- lasso.zeta
        file <- runjags::read.jagsfile(lasso.zeta)
        jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, a=a, b=b, alpha.precision=alpha.precision, coerce.precision=coerce.precision)
      }
      if (dim(w)[2]==0) {
        params <- c("alpha0", "beta", "lambda")
        inits[[1]]<-list("alpha0" = alpha, "beta" = rep(0,P))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "beta" = rep(0,P))
          }
        }
      } else {
        params <- c("alpha0", "beta", "lambda", "zeta")
        inits[[1]]<-list("alpha0" = alpha, "beta" = rep(0,P), "zeta" = rep(0, Q))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "beta" = rep(0,P), "zeta" = rep(0, Q))
          }
        }
      }
    } else if (model=="normalss")   {
      if (gamma.ind=="fixed") {
        if (dim(w)[2]==0) {
          file <- runjags::read.jagsfile(normalss.fixed)
          modelname <- normalss.fixed
          jagsdf <- list(y=y, X=x, N=N, P=P, k=k, pi.fixed=pi.fixed, alpha.precision=alpha.precision, sigma2.0=sigma2.0, sigma2.1=sigma2.1)
        } else {
          file <- runjags::read.jagsfile(normalss.fixed.zeta)
          modelname <- normalss.fixed.zeta
          jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, pi.fixed=pi.fixed, alpha.precision=alpha.precision, sigma2.0=sigma2.0, sigma2.1=sigma2.1, coerce.precision=coerce.precision)
        }
      } else {
        if (dim(w)[2]==0) {
          file <- runjags::read.jagsfile(normalss.random)
          modelname <- normalss.random
          jagsdf <- list(y=y, X=x, N=N, P=P, k=k, c.gamma=c.gamma, d.gamma=d.gamma, alpha.precision=alpha.precision, sigma2.0=sigma2.0, sigma2.1=sigma2.1)
        } else {
          file <- runjags::read.jagsfile(normalss.random.zeta)
          modelname <- normalss.random.zeta
          jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, c.gamma=c.gamma, d.gamma=d.gamma, alpha.precision=alpha.precision, sigma2.0=sigma2.0, sigma2.1=sigma2.1, coerce.precision=coerce.precision)
        }
      }
      if (dim(w)[2]==0) {
        params <- c("alpha0", "beta", "gamma")
        inits[[1]] <- list("alpha0" = alpha, "N1" = rep(0,P), "N2" = rep(0,P))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "N1" = rep(0,P), "N2" = rep(0,P))
          }
        }
      } else {
        params <- c("alpha0", "beta", "gamma", "zeta")
        inits[[1]] <- list("alpha0" = alpha, "N1" = rep(0,P), "N2" = rep(0,P), "zeta" = rep(0,Q))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "N1" = rep(0,P), "N2" = rep(0,P), "zeta" = rep(0,Q))
          }
        }
      }
    } else if (model=="dess") {
      if (gamma.ind=="fixed") {
        if (dim(w)[2]==0) {
          file <- runjags::read.jagsfile(dess.fixed)
          modelname <- dess.fixed
          jagsdf <- list(y=y, X=x, N=N, P=P, k=k, a=a, b=b, pi.fixed=pi.fixed, alpha.precision=alpha.precision, lambda0=lambda0)
        } else {
          file <- runjags::read.jagsfile(dess.fixed.zeta)
          modelname <- dess.fixed.zeta
          jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, a=a, b=b, pi.fixed=pi.fixed, alpha.precision=alpha.precision, lambda0=lambda0, coerce.precision=coerce.precision)
        }
      } else {
        if (dim(w)[2]==0) {
          file <- runjags::read.jagsfile(dess.random)
          modelname <- dess.random
          jagsdf <- list(y=y, X=x, N=N, P=P, k=k, a=a, b=b, c.gamma=c.gamma, d.gamma=d.gamma, alpha.precision=alpha.precision, lambda0=lambda0)
        } else {
          file <-runjags::read.jagsfile(dess.random.zeta)
          modelname <- dess.random.zeta
          jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, a=a, b=b, c.gamma=c.gamma, d.gamma=d.gamma, alpha.precision=alpha.precision, lambda0=lambda0, coerce.precision=coerce.precision)
        }
      }
      if (dim(w)[2]==0){
        params <- c("alpha0", "beta", "gamma", "lambda")
        inits[[1]] <- list("alpha0" = alpha, "DE1" = rep(0,P), "DE2" = rep(0,P))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "DE1" = rep(0,P), "DE2" = rep(0,P))
          }
        }
      } else {
        params <- c("alpha0", "beta", "gamma", "lambda", "zeta")
        inits[[1]] <- list("alpha0" = alpha, "DE1" = rep(0,P), "DE2" = rep(0,P), "zeta" = rep(0, Q))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "DE1" = rep(0,P), "DE2" = rep(0,P), "zeta" = rep(0, Q))
          }
        }
      }
    } else if (model=="regressvi") {
      if (gamma.ind=="fixed") {
        if (dim(w)[2]==0) {
          file <- runjags::read.jagsfile(regressvi.fixed)
          modelname <- regressvi.fixed
          jagsdf <- list(y=y, X=x, N=N, P=P, k=k, a=a, b=b, pi.fixed=pi.fixed, alpha.precision=alpha.precision)
        } else {
          file <- runjags::read.jagsfile(regressvi.fixed.zeta)
          modelname <- regressvi.fixed.zeta
          jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, a=a, b=b, pi.fixed=pi.fixed, alpha.precision=alpha.precision, coerce.precision=coerce.precision)
        }
      } else {
        if (dim(w)[2]==0) {
          file <- runjags::read.jagsfile(regressvi.random)
          modelname <- regressvi.random
          jagsdf <- list(y=y, X=x, N=N, P=P, k=k, a=a, b=b, c.gamma=c.gamma, d.gamma=d.gamma, alpha.precision=alpha.precision)

        } else {
          file <- runjags::read.jagsfile(regressvi.random.zeta)
          modelname <- regressvi.random.zeta
          jagsdf <- list(y=y, X=x, w=w, N=N, P=P, Q=Q, k=k, a=a, b=b, c.gamma=c.gamma, d.gamma=d.gamma, alpha.precision=alpha.precision, coerce.precision=coerce.precision)
        }
      }
      if (dim(w)[2]==0) {
        params <- c("alpha0", "gamma","betgma","lambda")
        inits[[1]] <-list("alpha0" = alpha, "beta" = rep(0,P))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "beta" = rep(0,P))
          }
        }
      } else {
        params <- c("alpha0", "gamma","betgma","lambda","zeta")
        inits[[1]] <-list("alpha0" = alpha, "beta" = rep(0,P), "zeta" = rep(0,Q))
        if (nChains > 1) {
          for (i in 2:nChains) {
            inits[[i]] <- list("alpha0" = sort(rnorm(k-1,0,.5)), "beta" = rep(0,P), "zeta" = rep(0,Q))
          }
        }
      }
    }
    if (parallel) {
          for (i in 1:length(inits)) {
            if (i%%2==0) {
              inits[[i]]$".RNG.name"<-'base::Wichmann-Hill'
              inits[[i]]$".RNG.seed"<-seed+i
            } else {
              inits[[i]]$".RNG.name"<-'base::Super-Duper'
              inits[[i]]$".RNG.seed"<-seed+i
            }
          }
    results<-runjags::run.jags(model=modelname, monitor=params, data = jagsdf, inits = inits,
          n.chains = nChains, adapt=adaptSteps, burnin=burnInSteps, sample=nIter,
          thin=thinSteps, method="parallel", n.sims=nChains)
    } else {
      if (is.null(seed)) seed<-26
      set.seed(seed)
      results<-runjags::run.jags(model=modelname, monitor=params, data = jagsdf, inits = inits,
                        n.chains = nChains, adapt=adaptSteps, burnin=burnInSteps, sample=nIter,
                        thin=thinSteps)
    }
    results<-runjags::add.summary(results)
    featureNames<-dimnames(x)[[2]]
    if (dim(w)[2]!=0) {
      #w<-as.matrix(object$w,drop=FALSE)
      varNames<-dimnames(w)[[2]]
      coda::varnames(results$mcmc)[grep("zeta",coda::varnames(results$mcmc))]<-varNames
      dimnames(results$summaries)[[1]][grep("zeta", dimnames(results$summaries)[[1]])]<-varNames
      dimnames(results$summary$statistics)[[1]][grep("zeta", dimnames(results$summary$statistics)[[1]])]<-varNames
      dimnames(results$summary$quantiles)[[1]][grep("zeta", dimnames(results$summary$statistics)[[1]])]<-varNames
      dimnames(results$HPD)[[1]][grep("zeta",dimnames(results$HPD)[[1]])]<-varNames
      dimnames(results$hpd)[[1]][grep("zeta",dimnames(results$hpd)[[1]])]<-varNames
 #     dimnames(results$mcse)[[1]][grep("zeta",dimnames(results$mcse)[[1]])]<-varNames
 #     colnames(results$autocorr)[grep("zeta",colnames(results$autocorr))]<-varNames
  #    dimnames(results$crosscorr)[[1]][grep("zeta",dimnames(results$crosscorr)[[1]])]<-varNames
 #     dimnames(results$crosscorr)[[2]][grep("zeta",dimnames(results$crosscorr)[[2]])]<-varNames
      names(results$truestochastic)[grep("zeta", names(results$truestochastic))]<-varNames
      names(results$semistochastic)[grep("zeta", names(results$semistochastic))]<-varNames
      names(results$nonstochastic)[grep("zeta", names(results$nonstochastic))]<-varNames
      names(results$discrete)[grep("zeta", names(results$discrete))]<-varNames
    }
    if (model=="regressvi" | model=="dess" | model=="normalss") {
      coda::varnames(results$mcmc)[grep("gamma",coda::varnames(results$mcmc))]<-paste("gamma", featureNames, sep=".")
      dimnames(results$summaries)[[1]][grep("gamma", dimnames(results$summaries)[[1]])]<-paste("gamma", featureNames, sep=".")
      dimnames(results$summary$statistics)[[1]][grep("gamma", dimnames(results$summary$statistics)[[1]])]<-paste("gamma", featureNames, sep=".")
      dimnames(results$summary$quantiles)[[1]][grep("gamma", dimnames(results$summary$statistics)[[1]])]<-paste("gamma", featureNames, sep=".")
      dimnames(results$HPD)[[1]][grep("gamma",dimnames(results$HPD)[[1]])]<-paste("gamma", featureNames, sep=".")
      dimnames(results$hpd)[[1]][grep("gamma",dimnames(results$hpd)[[1]])]<-paste("gamma", featureNames, sep=".")
 #     dimnames(results$mcse)[[1]][grep("gamma",dimnames(results$mcse)[[1]])]<-paste("gamma", featureNames, sep=".")
 #     fnumber <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(results$autocorr)))
 #     colnames(results$autocorr)[grep("gamma",colnames(results$autocorr))]<-paste("gamma", featureNames[fnumber[grep("bet",colnames(results$autocorr))]], sep=".")
 #     dimnames(results$crosscorr)[[1]][grep("gamma",dimnames(results$crosscorr)[[1]])]<-paste("gamma", featureNames, sep=".")
 #     dimnames(results$crosscorr)[[2]][grep("gamma",dimnames(results$crosscorr)[[2]])]<-paste("gamma", featureNames, sep=".")
      names(results$truestochastic)[grep("gamma",names(results$truestochastic))]<-paste("gamma", featureNames, sep=".")
      names(results$semistochastic)[grep("gamma",names(results$semistochastic))]<-paste("gamma", featureNames, sep=".")
      names(results$nonstochastic)[grep("gamma",names(results$nonstochastic))]<-paste("gamma", featureNames, sep=".")
      names(results$discrete)[grep("gamma",names(results$discrete))]<-paste("gamma", featureNames, sep=".")
    }
    coda::varnames(results$mcmc)[grep("bet",coda::varnames(results$mcmc))]<-paste("beta", featureNames, sep=".")
    dimnames(results$summaries)[[1]][grep("bet", dimnames(results$summaries)[[1]])]<-paste("beta", featureNames, sep=".")
    dimnames(results$summary$statistics)[[1]][grep("bet", dimnames(results$summary$statistics)[[1]])]<-paste("beta", featureNames, sep=".")
    dimnames(results$summary$quantiles)[[1]][grep("bet", dimnames(results$summary$statistics)[[1]])]<-paste("beta", featureNames, sep=".")
    dimnames(results$HPD)[[1]][grep("bet",dimnames(results$HPD)[[1]])]<-paste("beta", featureNames, sep=".")
    dimnames(results$hpd)[[1]][grep("bet",dimnames(results$hpd)[[1]])]<-paste("beta", featureNames, sep=".")
 #   dimnames(results$mcse)[[1]][grep("bet",dimnames(results$mcse)[[1]])]<-paste("beta", featureNames, sep=".")
 #   fnumber <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(results$autocorr)))
  #  colnames(results$autocorr)[grep("bet",colnames(results$autocorr))]<-paste("beta", featureNames[fnumber[grep("bet",colnames(results$autocorr))]], sep=".")
    #dimnames(results$crosscorr)[[1]][grep("bet",dimnames(results$crosscorr)[[1]])]<-paste("beta", featureNames, sep=".")
    #dimnames(results$crosscorr)[[2]][grep("bet",dimnames(results$crosscorr)[[2]])]<-paste("beta", featureNames, sep=".")
    names(results$truestochastic)[grep("bet", names(results$truestochastic))]<-paste("beta", featureNames, sep=".")
    names(results$semistochastic)[grep("bet", names(results$semistochastic))]<-paste("beta", featureNames, sep=".")
    names(results$nonstochastic)[grep("bet", names(results$nonstochastic))]<-paste("beta", featureNames, sep=".")
    names(results$discrete)[grep("bet", names(results$discrete))]<-paste("beta", featureNames, sep=".")
    if (model=="lasso") {
      output<-list(results, cl, model, a, b, featureNames, center, scale, Y, x, w)
      names(output)<-c("results", "call", "model", "a", "b","featureNames", "center", "scale", "y", "x", "w")
    } else if (model=="regressvi") {
      output<-list(results, cl, model, a, b, gamma.ind, pi.fixed, c.gamma, d.gamma, featureNames, center, scale, Y, x, w)
      names(output)<-c("results", "call", "model", "a", "b", "gamma.ind", "pi.fixed", "c.gamma", "d.gamma", "featureNames", "center", "scale", "y", "x", "w")
    } else if (model=="normalss") {
      output<-list(results, cl, model, a, b, gamma.ind, pi.fixed, c.gamma, d.gamma, sigma2.0, sigma2.1, featureNames, center, scale, Y, x, w)
      names(output)<-c("results", "call", "model", "a", "b", "gamma.ind", "pi.fixed", "c.gamma", "d.gamma", "sigma2.0", "sigma2.1", "featureNames", "center", "scale", "y", "x", "w")
    } else if (model=="dess") {
      output<-list(results, cl, model, a, b, gamma.ind, pi.fixed, c.gamma, d.gamma, lambda0, featureNames, center, scale, Y, x, w)
      names(output)<-c("results", "call", "model", "a", "b", "gamma.ind", "pi.fixed", "c.gamma", "d.gamma", "lambda0", "featureNames", "center", "scale", "y", "x", "w")
    }
    class(output)<-"ordinalbayes"
    output
  }
