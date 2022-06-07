library(fastGHQuad)

#######Functions used to maximize missing data likelihood#######
### Likelihood Function ###
Like <- function(eta, Y, Cov, sigma, alpha, delta, t.df.miss=4, rule) {
  obs <- !is.na(Y)
  mu <- as.vector(Cov%*%eta) 
  g1 <- sapply(X = which(!obs), function(i) {
    fastGHQuad::ghQuad(f = g1.quad, rule = rule, mu=mu[i], sigma=sigma, delta=delta, alpha=alpha, t.df.miss=t.df.miss)
  })
  return(-1/2/sigma^2*sum((Y[obs]-mu[obs])^2) + sum(log(1-g1)))
}

g1.quad <- function(x, mu, sigma, delta, alpha, t.df.miss=4) {
  1/sqrt(pi) * pt(q = alpha*(mu+sqrt(2)*sigma*x-delta), df=t.df.miss)
}

## Gradient ##
Grad <- function(eta, Y, Cov, sigma, alpha, delta, t.df.miss=4, rule) {
  obs <- !is.na(Y)
  mu <- as.vector(Cov%*%eta)
  g1 <- sapply(X = which(!obs), function(i) {
    fastGHQuad::ghQuad(f = g1.quad, rule = rule, mu=mu[i], sigma=sigma, delta=delta, alpha=alpha, t.df.miss=t.df.miss)
  })
  g2 <- sapply(X = which(!obs), function(i) {
    fastGHQuad::ghQuad(f = g2.quad, rule = rule, mu=mu[i], sigma=sigma, delta=delta, alpha=alpha, t.df.miss=t.df.miss)
  })
  s <- rep(NA, length(mu)) 
  s[ obs] <- 1/sigma^2*(Y[obs]-mu[obs]) 
  s[!obs] <- -alpha*g2/(1-g1) 
  return(as.vector(t(Cov)%*%s)) 
}


g2.quad <- function(x, mu, sigma, delta, alpha, t.df.miss=4) {
  1/sqrt(pi) * dt(x = alpha*(mu+sqrt(2)*sigma*x-delta), df=t.df.miss)
}

#To be used to estimate the variance
Covn <- function(eta, Y, Cov, sigma, alpha, delta, t.df.miss=4, rule, use.empirical=F, return.all=F) {
  obs <- !is.na(Y)
  mu <- as.vector(Cov%*%eta)
  n <- length(mu)
  g1 <- sapply(X = mu, function(mu.i) {
    fastGHQuad::ghQuad(f = g1.quad, rule = rule, mu=mu.i, sigma=sigma, delta=delta, alpha=alpha, t.df.miss=t.df.miss)
  })
  g2 <- sapply(X = mu, function(mu.i) {
    fastGHQuad::ghQuad(f = g2.quad, rule = rule, mu=mu.i, sigma=sigma, delta=delta, alpha=alpha, t.df.miss=t.df.miss)
  })
  g3 <- sapply(X = mu, function(mu.i) {
    fastGHQuad::ghQuad(f = g3.quad, rule = rule, mu=mu.i, sigma=sigma, delta=delta, alpha=alpha, t.df.miss=t.df.miss)
  })
  ai <- g1/sigma^2 + alpha^2*(g3+g2^2/(1-g1))
  s <- -alpha*g2/(1-g1)
  s[obs] <- 1/sigma^2*(Y[obs]-mu[obs])
  hessian <- (-t(Cov*ai)%*%Cov)
  Hess.inverse <- solve(hessian)
  if (return.all) {out <- list(); out$NaiveCov <- -Hess.inverse}
  if (!use.empirical) {
    v <- 1/( 1 + ai*rowSums((Cov%*%Hess.inverse)*Cov) )^2
    if (return.all) {out$Cov <- Hess.inverse %*% t(Cov*s*v)%*%(Cov*s) %*% Hess.inverse; return(out)}
    return(Hess.inverse %*% t(Cov*s*v)%*%(Cov*s) %*% Hess.inverse)
  } else {
    ai.hat <- rep(1/sigma^2,length(obs)); ai.hat[!obs] <- pmax(0,alpha^2*(g3+g2^2/(1-g1))/(1-g1))[!obs]
    v <- 1/( 1 - ai.hat*rowSums((Cov%*%solve( t(Cov*ai.hat)%*%Cov ))*Cov) )^2
    if (return.all) {out$Cov <- Hess.inverse %*% t(Cov*s*v)%*%(Cov*s) %*% Hess.inverse; return(out)}
    return(Hess.inverse %*% t(Cov*s*v)%*%(Cov*s) %*% Hess.inverse)
  }
}

g3.quad <- function(x, mu, sigma, delta, alpha, t.df.miss=4) {
  x <- alpha*(mu + x*sqrt(2)*sigma - delta)
  1/sqrt(pi) * dt(x = x, df=t.df.miss)*(-x*(t.df.miss+1)/(t.df.miss+x^2))
}
