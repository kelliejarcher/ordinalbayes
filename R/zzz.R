regressvi.fixed <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],betgma[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (m in 1:P){
      beta[m] ~ ddexp(0,lambda)
    }
    lambda ~ dgamma(a, b)
    for (j in 1:P){
      betgma[j] <- beta[j] * gamma[j]
      gamma[j] ~ dbern(pi.fixed)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

regressvi.fixed.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],betgma[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(l in 2:(k-1)){
        logit(U[i,l]) <- alpha[l]-mu[i]
        p[i,l] <- U[i,l] - U[i,l-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (m in 1:P){
      beta[m] ~ ddexp(0,lambda)
    }
    lambda ~ dgamma(a, b)
    for (j in 1:P){
      betgma[j] <- beta[j] * gamma[j]
      gamma[j] ~ dbern(pi.fixed)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

regressvi.random <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],betgma[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (m in 1:P){
      beta[m] ~ ddexp(0,lambda)
    }
    lambda ~ dgamma(a, b)
    for (j in 1:P){
      betgma[j] <- beta[j] * gamma[j]
      gamma[j] ~ dbern(pi[j])
      pi[j] ~ dbeta(c.gamma, d.gamma)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

regressvi.random.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],betgma[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(l in 2:(k-1)){
        logit(U[i,l]) <- alpha[l]-mu[i]
        p[i,l] <- U[i,l] - U[i,l-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (m in 1:P){
      beta[m] ~ ddexp(0,lambda)
    }
    lambda ~ dgamma(a, b)
    for (j in 1:P){
      betgma[j] <- beta[j] * gamma[j]
      gamma[j] ~ dbern(pi[j])
      pi[j] ~ dbeta(c.gamma, d.gamma)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

dess.fixed <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],beta[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (j in 1:P){
      beta[j] <- DE1[j] * (1 - gamma[j]) + DE2[j] * gamma[j]
      gamma[j] ~ dbern(pi.fixed)
      DE1[j] ~ ddexp(0, lambda0) # variance = 0.005
      DE2[j] ~ ddexp(0, lambda) # lambda is assigned a gamma prior
    }
    lambda ~ dgamma(a, b)
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

dess.fixed.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],beta[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(l in 2:(k-1)){
        logit(U[i,l]) <- alpha[l]-mu[i]
        p[i,l] <- U[i,l] - U[i,l-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (j in 1:P){
      beta[j] <- DE1[j] * (1 - gamma[j]) + DE2[j] * gamma[j]
      DE1[j] ~ ddexp(0, lambda0)
      DE2[j] ~ ddexp(0, lambda)
      gamma[j] ~ dbern(pi.fixed)
    }
    lambda ~ dgamma(a, b)
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

dess.random <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],beta[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (j in 1:P){
      beta[j] <- DE1[j] * (1 - gamma[j]) + DE2[j] * gamma[j]
      gamma[j] ~ dbern(pi[j])
      pi[j] ~ dbeta(c.gamma, d.gamma)
      DE1[j] ~ ddexp(0, lambda0) # variance = 0.005
      DE2[j] ~ ddexp(0, lambda) # lambda is assigned a gamma prior
    }
    lambda ~ dgamma(a, b)
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

dess.random.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],beta[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(l in 2:(k-1)){
        logit(U[i,l]) <- alpha[l]-mu[i]
        p[i,l] <- U[i,l] - U[i,l-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (j in 1:P){
      beta[j] <- DE1[j] * (1 - gamma[j]) + DE2[j] * gamma[j]
      DE1[j] ~ ddexp(0, lambda0)
      DE2[j] ~ ddexp(0, lambda)
      gamma[j] ~ dbern(pi[j])
      pi[j] ~ dbeta(c.gamma, d.gamma)
    }
    lambda ~ dgamma(a, b)
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

lasso <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],beta[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (m in 1:P){
      beta[m] ~ ddexp(0,lambda)
    }
    lambda ~ dgamma(a, b)
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

lasso.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],beta[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(j in 2:(k-1)){
        logit(U[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- U[i,j] - U[i,j-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (m in 1:P){
      beta[m] ~ ddexp(0,lambda)
    }
    lambda ~ dgamma(a, b)
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

normalss.fixed <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],beta[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (j in 1:P){
      beta[j] <- N1[j] * (1 - gamma[j]) + N2[j] * gamma[j]
      gamma[j] ~ dbern(pi.fixed)
      N1[j] ~ dnorm(0, 1/sigma2.0)
      N2[j] ~ dnorm(0, 1/sigma2.1)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

normalss.fixed.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],beta[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(l in 2:(k-1)){
        logit(U[i,l]) <- alpha[l]-mu[i]
        p[i,l] <- U[i,l] - U[i,l-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (j in 1:P){
      beta[j] <- N1[j] * (1 - gamma[j]) + N2[j] * gamma[j]
      N1[j] ~ dnorm(0, 1/sigma2.0)
      N2[j] ~ dnorm(0, 1/sigma2.1)
      gamma[j] ~ dbern(pi.fixed)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

normalss.random <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(X[i,],beta[])
      logit(Q[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- Q[i,1]
      for(j in 2:(k-1)){
        logit(Q[i,j]) <- alpha[j]-mu[i]
        p[i,j] <- Q[i,j] - Q[i,j-1]
      }
      p[i,k] <- 1 - Q[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (j in 1:P){
      beta[j] <- N1[j] * (1 - gamma[j]) + N2[j] * gamma[j]
      gamma[j] ~ dbern(pi[j])
      pi[j] ~ dbeta(c.gamma, d.gamma)
      N1[j] ~ dnorm(0, 1/sigma2.0)
      N2[j] ~ dnorm(0, 1/sigma2.1)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"

normalss.random.zeta <- "
  model{
    for(i in 1:N){
      mu[i] <- inprod(w[i,],zeta[])+inprod(X[i,],beta[])
      logit(U[i,1]) <- alpha[1]-mu[i]
      p[i,1] <- U[i,1]
      for(l in 2:(k-1)){
        logit(U[i,l]) <- alpha[l]-mu[i]
        p[i,l] <- U[i,l] - U[i,l-1]
      }
      p[i,k] <- 1 - U[i,(k-1)]
      y[i,] ~ dcat(p[i,1:k])
    }
    for (l in 1:Q){
      zeta[l] ~ dnorm(0, coerce.precision)
    }
    for (j in 1:P){
      beta[j] <- N1[j] * (1 - gamma[j]) + N2[j] * gamma[j]
      N1[j] ~ dnorm(0, 1/sigma2.0)
      N2[j] ~ dnorm(0, 1/sigma2.1)
      gamma[j] ~ dbern(pi[j])
      pi[j] ~ dbeta(c.gamma, d.gamma)
    }
    for (r in 1:(k-1)) {
      alpha0[r] ~ dnorm(0, alpha.precision)
    }
    alpha <- sort(alpha0)
  }
"



