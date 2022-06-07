library(bcv)
library(impute)
library(missForest)
library(cate)

### Impute dataset with missing values ###
## data: #metabolites by #samples data matrix
## methods: a vector containing methods you selected
## max.miss: the maximum missing proportion of metabolites that will be included in analysis

impute.missing <- function (
  data,
  methods = c("MIN","SVD","KNN","RF"),
  max.miss = 0.5
) {
  
  p <- nrow(data)
  n <- ncol(data)
  
  Frac.missing <- rowMeans(is.na(data))
  ind.mnar <- which(Frac.missing <= max.miss)
  
  results <- list()
  
  ## min value among variable
  if("MIN" %in% methods)  {
    cat("Imputing by MIN...")
    out <- data
    for (j in ind.mnar) out[j,][is.na(out[j,])] <- min(out[j,],na.rm=T)
    results <- c(results, MIN=list(out))
    cat("done\n")}
  
  ## SVD decomposition
  if ("SVD" %in% methods)  {
    cat("Imputing by SVD...")
    out <- data
    out[ind.mnar,] <- t(bcv::impute.svd(x = t(data[ind.mnar,]), k = 10)$x)
    results <- c(results, SVD=list(out))
    cat("done\n")}
  
  # KNN on metabolites
  if ("KNN" %in% methods) {
    cat("Imputing by KNN...")
    out <- data
    out[ind.mnar,] <- t(impute::impute.knn(data = t(data[ind.mnar,]))$data)
    results <- c(results, KNN=list(out))
    cat("done\n")}
  
  if ("RF" %in% methods) {
    out <- data
    out[ind.mnar,] <- t(missForest::missForest(xmis = t(data[ind.mnar,]), maxiter = 10)$ximp)
    results <- c(results, RF=list(out))
  }
  
  return(results)
}


### Get latent factors estimation and OLS estimation ###
## X: the phenotype of interest.
## Z: the observed nuisance covariates.
## Y: the #metabolites x #samples imputed data matrix.
## K: the number of latent factors. Should be the same for all methods.

get.OLS <- function(Y, X, Z, K, missing.max=0.5) {
  if(!is.null(Z)) Z <- cbind(Z)
  if(!is.null(X)) X <- cbind(X)
  
  p <- ifelse(is.null(Y),0,NROW(Y))
  r <- ifelse(is.null(Z),0,NCOL(Z))
  d <- ifelse(is.null(X),0,NCOL(X))
  
  Frac.Missing <- rowMeans(is.na(Y))
  ind.obs <- which(Frac.Missing == 0)
  
  #### Function to get latent variable using SVD ####
  Latent.SVD <- function(Y, K) {
    n <- ncol(Y)
    Z <- cbind(rep(1, n))
    P <- diag(n) - Z %*% solve(t(Z)%*%Z)%*%t(Z)
    tmp <- svd(Y %*% P, nu = K, nv = K)
    out <- list()
    out$C <- tmp$v %*% diag(sqrt(tmp$d[1:K]))/sqrt(n-1)
    out$L <- tmp$u %*% diag(tmp$d[1:K]) / sqrt(n-1)
    return(out)
  }
  
  # Estimate C
  C <- cate::cate.fit(X.primary = cbind(X), X.nuis =cbind(Z), Y = t(Y[ind.obs,]), r = K, calibrate = F)$Z
  
  Cov <- cbind(X,Z,C)
  
  B <- matrix(data = NA, nrow = p, ncol = d)
  L <- matrix(data = NA, nrow = p, ncol = K)
  Sigma <- rep(NA, p)
  Var.Beta <- matrix(data = NA, nrow = p, ncol = d)
  Var.list <- vector(mode = "list", length = nrow(Y))
  p.values <- matrix(data = NA, nrow = p, ncol = d)
  
  # OLS on metabolites with no missing value
  B.all <- Y[ind.obs,]%*%Cov%*%solve(t(Cov)%*%Cov)
  B[ind.obs,] <- B.all[,1:d]
  L[ind.obs,] <- B.all[,(d+r+1):(d+r+K)]
  Sigma[ind.obs] <- rowSums((Y[ind.obs,]-B.all%*%t(Cov))^2)/(nrow(X)-ncol(Cov))
  Var.Beta[ind.obs,] <- cbind(Sigma[ind.obs]) %*% diag(solve(t(Cov)%*%Cov))[1:d]
  Var.list[ind.obs ] <- lapply(ind.obs, function(i){Sigma[i] * solve(t(Cov)%*%Cov)})
  p.values[ind.obs,] <- 2*pnorm(-abs( B[ind.obs,] / sqrt(Var.Beta[ind.obs,]) ))
  
  # OLS on metabolites with missing values
  ind.miss <- which(Frac.Missing > 0 & Frac.Missing <= missing.max)
  if (length(ind.miss) > 0) {
    for (g in ind.miss) {
      ind.g <- !is.na(Y[g,])
      y <- cbind(Y[g,ind.g])
      x <- cbind(Cov[ind.g,])
      hess <- solve(t(x) %*% x)
      beta.hat <- hess %*% t(x) %*% y
      residual <- as.vector(y - x %*% beta.hat)
      B[g,] <- beta.hat[1:d]
      L[g,] <- beta.hat[-(1:(d+r))]
      Sigma[g] <- sum(residual^2) / (length(y) - NCOL(x))
      Var.Beta[g,] <- Sigma[g]*diag(hess)[1:d]
      p.values[g,] <- 2*pnorm(-abs( B[g,] / sqrt(Var.Beta[g,]) ))
    }
  }
  return(list(B=B, C=C, L=L, Sigma=Sigma, Var.Beta=Var.Beta, Var.list=Var.list, p.t=p.values))
}

### Compute Confidence Interval Coverage ###
## Z: Z scores of estimated parameters of interest
## B: simulated parameters of interest
## Frac: missing fraction of metabolites

get.CIC <- function(Z, B, Frac, breaks, min.miss=0.05, max.miss=0.5, alpha=0.05) {
  CI.coverage <- rep(NA, length(breaks)-1)
  for (j in 1:(length(B.breaks)-1)) {
    ind <- abs(B) > B.breaks[j] & abs(B) <= B.breaks[j+1] & Frac > min.miss & Frac <= max.miss
    tmp <- Z[ind] >= qnorm(alpha/2) & Z[ind] <= qnorm(1-alpha/2)
    CI.coverage[j] <- mean(tmp, na.rm = T)
  }
  return(CI.coverage)
}

### Compute Confidence Interval Width ###
## V: estimated variance of parameters of interest
## B: simulated parameters of interest
## Frac: missing fraction of metabolites

get.CIW <- function(V, B, Frac, breaks, min.miss=0.05, max.miss=0.5, alpha=0.05) {
  out <- list()
  out$breaks <- c()
  out$CIW <- c()
  for (j in 1:(length(B.breaks)-1)) {
    ind <- abs(B) > B.breaks[j] & abs(B) <= B.breaks[j+1] & Frac > min.miss & Frac <= max.miss
    out$breaks <- c(out$breaks, rep(j, sum(ind)))
    ciw <- na.omit(2*qnorm(1-alpha/2) * sqrt(V[ind]))
    out$CIW <- c(out$CIW, ciw)
  }
  return(out)
}

### FDP ###

sim.FDP <- function(p.values, Beta, Frac.Missing, Ind.Confident = NULL, ii) {
  FDP <- matrix(NA, ii, 4)
  for (i in 1:ii) {
    idx.all <- Frac.Missing[,i] <= 0.05 | Ind.Confident[,i]
    FDP[i,] <- get.FDP(p.values = p.values[idx.all,i], Beta = Beta[idx.all,i])
  }
  return(FDP)
}

get.FDP <- function(p.values, Beta, q.points = seq(0.05, 0.2, by=0.05), pi0=NULL) {
  Beta <- Beta[!is.na(p.values)]
  p.values <- p.values[!is.na(p.values)]
  if (is.null(pi0)) {
    q.values <- (qvalue::qvalue(p.values))$qvalue
  } else {
    q.values <- (qvalue::qvalue(p.values, pi0 = pi0))$qvalue
  }
  
  FDP <- rep(NA, length(q.points))
  for (i in 1:length(q.points)) {
    FDP[i] <- sum(abs(Beta) < 1e-8 & q.values <= q.points[i]) / sum(q.values <= q.points[i])
  }
  return(FDP)
}

### Power ###

sim.Power <- function(p.values, Beta, Frac.Missing, Ind.Confident) {
  Power <- matrix(NA, ncol(Beta), 4)
  for (i in 1:n.sim) {
    idx.all <- Frac.Missing[,i] <= 0.05 | Ind.Confident[,i]
    Power[i,] <- get.Power(p.values = p.values[idx.all,i], Beta = Beta[idx.all,i])
  }
  return(Power)
}

get.Power <- function(p.values, Beta, idx.mis, q.points = seq(0.05, 0.20, by=0.05), pi0=NULL) {
  Beta <- Beta[!is.na(p.values)]
  p.values <- p.values[!is.na(p.values)]
  if (is.null(pi0)) {
    q.values <- (qvalue::qvalue(p.values))$qvalue
  } else {
    q.values <- (qvalue::qvalue(p.values, pi0 = pi0))$qvalue
  }
  
  Power <- rep(NA, length(q.points))
  for (i in 1:length(q.points)) {
    Power[i] <- sum(abs(Beta) >= 1e-8 & q.values <= q.points[i]) / sum(abs(Beta) >= 1e-8)
  }
  return(Power)
}