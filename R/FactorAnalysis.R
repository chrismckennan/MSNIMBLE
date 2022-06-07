library(fastGHQuad)
library(sva)

#' Factor analysis with non-random missing data
#' 
#' Perform factor analysis in metabolite data with non-random missing observations. The goal is to estimate L and C in the model Y = LC' + E, where we identify L and C using "PCA identification", i.e. C'C is proportional to the identity and L'L is diagonal with decreasing elements.
#' 
#' @param Y a \code{p} x \code{n} data matrix of log2-transformed metabolite intensities, where \code{p} = #of metabolites and \code{n} = #of samples. Missing values should be left as \code{NA}.
#' @param Cov a \code{n} x \code{r} matrix of nuisance covariates (i.e. intercept, sex, etc.).
#' @param K The number of latent covariates (i.e. C is a \code{n} x \code{K} matrix). If unspecified, it is estimated using sva::num.sv applied to the metabolites with complete data.
#' @param Miss.Mech The missingness mechansim object returned by \code{EstimateMissing}.
#' @param ind.samples A logical or numeric vector of samples to be considered in the analysis. For example, if nuisance covariates were only measured in a subset of the samples, this would be the samples with nuisance covariates. Default is no missing samples.
#' @param method The method used to estimate L. Can be one of "MSNIMBLE" (the default and recommended; uses both observed and missing data to improve power), "sIPW" (stabilized inverse probability weighting; discards missing data), or "IPW" (inverse probability weighting; discards missing data).
#' @param max.miss Maximum fraction of missing data a metabolite is allowed to have to be used in the estimate for the latent factors C. Defaults to 0.5.
#' @param n.quad Number of Gaussian quadrature points. Default, and recommended value, is 100.
#' 
#' @return A list \item{C}{Estimates for the latent factors. The columns are orthogonal, orthogonal to Cov, and have squared norm \code{NROW(Cov)-NCOL(Cov)}.} \item{L}{Estimates for the latent loadings. Columns are orthogonal, and are arranged so that their squared norms are decreasing.} \item{Sigma}{The metabolite-specific residuals variances.} \item{Var.L}{A list of K x K variances of the rows of L.} \item{stand.eigen.values}{Standardized eigenvalue estimates of 1/#metabolites*L(C'C)L', which have been standardized \code{mean(Sigma)}.} \item{Perc.Var.Explained}{The percent metabolite level variance explained by the \code{K} latent factors.} \item{method}{The method used to estimate L.}
#' @export
FactorAnalysis <- function(Y, Cov=NULL, K=NULL, Miss.Mech, ind.samples=NULL, method=c("MSNIMBLE", "sIPW", "IPW"), max.miss=0.5, BH.min=NULL, include.updates=T, n.quad=100) {
  out <- list()
  method <- match.arg(method,c("MSNIMBLE", "sIPW", "IPW"))
  Alpha <- Miss.Mech$Post.Theta[,1]
  Delta <- Miss.Mech$Post.Theta[,2]
  max.miss.C <- Miss.Mech$max.miss.C
  t.df.miss <- ifelse(is.null(Miss.Mech$t.df.miss),4,Miss.Mech$t.df.miss)
  if (is.null(Cov)) {Cov <- rep(1,ncol(Y))}
  Cov <- cbind(Cov)
  
  #Fix observed samples#
  if (is.null(ind.samples)) {
    ind.samples <- 1:ncol(Y)
  } else {
    if (is.logical(ind.samples[1])) {
      length.samples <- sum(ind.samples)
    } else {
      length.samples <- length(ind.samples)
    }
    if (ncol(Y) > length.samples) {Y <- Y[,ind.samples]}
    if (nrow(Cov) > length.samples) {Cov <- cbind(Cov[ind.samples,])}
  }
  
  #Preliminary objects#
  d <- ncol(Cov)
  p <- nrow(Y)
  n <- ncol(Y)
  Prob.Missing <- rowMeans(is.na(Y))
  out$Cov <- Cov; out$method <- method
  ind.miss.all <- !is.na(Miss.Mech$Theta.Miss[,1]) & !is.na(Miss.Mech$Theta.Miss[,2])  #Metabolites with missing data
  if (is.null(BH.min)) {BH.min <- Miss.Mech$BH.min}
  tmp <- BH.proc(p = Miss.Mech$Pvalue.value, alpha = BH.min)
  
  out$flag <- !(tmp == FALSE & !is.na(tmp)) & Prob.Missing > max.miss.C  #All flagged metabolites
  out$flag.analyzed <- !(tmp == FALSE & !is.na(tmp)) & ind.miss.all    #Metabolites that have reported p-values, but the missingness mechanism is suspect
  
  
  #Determine Weights#
  if (method%in%c("sIPW","MSNIMBLE")) {
    Weights <- Miss.Mech$Post.W[,ind.samples] * Miss.Mech$Pi.MAR[,ind.samples]
    Var.Weights <- Miss.Mech$Pi.MAR[,ind.samples]^2*( Miss.Mech$Post.W[,ind.samples]^2 + Miss.Mech$Post.VarW[,ind.samples] )
  }
  if (method == "IPW") {
    Weights <- Miss.Mech$Post.W[,ind.samples]
    Var.Weights <- Miss.Mech$Post.W[,ind.samples]^2 + Miss.Mech$Post.VarW[,ind.samples]
  }
  Weights[is.na(Weights)] <- 1
  
  #Get K#
  if (is.null(K)) {
    if (include.updates) {cat("Estimating K ...")}
    K <- sva::num.sv(dat = Y[Prob.Missing==0,], mod = cbind(out$X,out$Z))
    if (include.updates) {cat("done\n")}
  }; out$K <- K
  
  #Preliminary estimate for C#
  if (include.updates) {cat(paste0("Estimating ", K, " latent factors..."))}
  out.C <- EstC.0(Y = Y, K = K, Cov = cbind(Cov), max.miss = max.miss.C, max.iter = 800, n.repeat.Sigma = 1)
  out$C <- out.C$C
  if (include.updates) {cat("done\n")}
  
  #Refine estimate for C.perp#
  if (max.miss > max.miss.C) {
    if (include.updates) {cat(paste("Refining estimate for latent factors with missingness mech...", collapse = ""))}
    ind.perp <- Prob.Missing <= max.miss.C | (Miss.Mech$Ind.Confident & Prob.Missing <= max.miss)
    out.perp <- Est.Cperp.Weights(Y = Y, Cov = cbind(Cov), C.start = out$C, Weights = Weights, ind.use.miss = ind.perp, max.miss.C = max.miss.C)
    out$n.iter.refine <- out.perp$n.iter
    out$C <- out.perp$C.perp
    if (include.updates) {cat("done\n")}
  } else {
    ind.perp <- Prob.Missing <= max.miss.C
  }
  
  #Get estimates for L#
  ind.obs <- which(Prob.Missing == 0)
  ind.mar <- which(Prob.Missing > 0 & Prob.Missing <= max.miss.C)
  ind.mnar <- which(Prob.Missing > max.miss.C & ind.miss.all)
  Cov.total <- cbind(out$C,Cov)
  LB <- matrix(NA,nrow=p,ncol=K+d)
  out$Sigma <- rep(NA,p)
  out$Var.L <- vector(mode = "list", length = p)
  if (length(ind.obs)==1) {ind.mar <- c(ind.mar,ind.obs)}
  if (include.updates) {cat("Estimating loadings for metabolites with complete data...")}
  if (length(ind.obs) > 1) {
    LB[ind.obs,] <- Y[ind.obs,]%*%Cov.total%*%solve(t(Cov.total)%*%Cov.total)
    out$Sigma[ind.obs] <- rowSums((Y[ind.obs,] - LB[ind.obs,]%*%t(Cov.total))^2)/(n-K-d)
    tmp <- solve(t(Cov.total)%*%Cov.total)[1:K,1:K]
    out$Var.L[ind.obs] <- lapply(ind.obs,function(g){out$Sigma[g]*tmp}); rm(tmp)
  }
  if (length(ind.mar) > 0) {
    tmp <- lapply(ind.mar,function(g){
      ind.obs <- !is.na(Y[g,])
      y <- Y[g,ind.obs]
      cov.tmp <- Cov.total[ind.obs,]
      hess <- solve(t(cov.tmp)%*%cov.tmp)
      tt <- list(); tt$coef <- hess%*%t(cov.tmp)%*%y
      tt$sigma2 <- 1/(length(y)-ncol(cov.tmp))*sum((y-cov.tmp%*%tt$coef)^2)
      tt$var.C <- tt$sigma2*hess[1:K,1:K]
      return(tt)
    })
    LB[ind.mar,] <- t(sapply(tmp,function(x){x$coef}))
    out$Sigma[ind.mar] <- sapply(tmp,function(x){x$sigma2})
    out$Var.L[ind.mar] <- lapply(tmp,function(x){x$var.C}); rm(tmp)
  }
  if (include.updates) {cat("done\n")}
  if (include.updates) {cat("Estimating loadings for metabolites with missing data...")}
  if (length(ind.mnar) > 0) {
    if (method=="MSNIMBLE") {rule.gh <- fastGHQuad::gaussHermiteData(n.quad)}
    tmp <- lapply(ind.mnar,function(g){
      ind.obs <- !is.na(Y[g,])
      y <- Y[g,ind.obs]
      w <- Weights[g,ind.obs]
      var.w <- Var.Weights[g,ind.obs]
      cov.tmp <- cbind(Cov.total[ind.obs,])
      hess <- solve(t(cov.tmp*w)%*%cov.tmp)
      tt <- list(); tt$coef <- hess%*%t(cov.tmp*w)%*%y
      resids.ipw <- y-as.vector(cov.tmp%*%tt$coef)
      hat.ipw <- apply(X = cov.tmp*sqrt(w), MARGIN = 1, function(x,A){sum(x*as.vector(A%*%x))}, A=hess)
      tt$sigma2 <- sum( w*resids.ipw^2 / (1-hat.ipw)^2 )/sum(w)
      if (method=="MSNIMBLE") {
        out.g <- optim(par=tt$coef, fn=Like, gr=Grad, Y=Y[g,], Cov=Cov.total, sigma=sqrt(tt$sigma2), alpha=Alpha[g], delta=Delta[g], t.df.miss=t.df.miss, rule=rule.gh, control=list(fnscale=-1), method="BFGS")
        tt$coef <- out.g$par
        out.cov.g <- Covn(eta = tt$coef, Y = Y[g,], Cov = Cov.total, sigma = sqrt(tt$sigma2), alpha = Alpha[g], delta = Delta[g], t.df.miss = t.df.miss, rule = rule.gh, use.empirical = F, return.all = T)
        tt$var.vc <- out.cov.g$Cov
        return(tt)
      }
      Score <- cov.tmp*resids.ipw
      Score.vc <- t(sapply(X = 1:nrow(cov.tmp), function(i){ as.vector(solve( diag(1,ncol(cov.tmp),ncol(cov.tmp)) - w[i]*cbind(cov.tmp[i,])%*%rbind(cov.tmp[i,])%*%hess, Score[i,] )) }))
      tt$var.vc <- hess%*%( t(Score.vc*var.w)%*%Score.vc )%*%hess
      return(tt)
    })
    LB[ind.mnar,] <- t(sapply(tmp,function(x){x$coef}))
    out$Sigma[ind.mnar] <- sapply(tmp,function(x){x$sigma2})
    out$Var.L[ind.mnar] <- lapply(tmp,function(x){x$var.vc[1:K,1:K]}); rm(tmp)
  }
  if (include.updates) {cat("done\n")}
  
  out$L <- LB[,1:K]
  R <- chol(1/(n-d)*t(out$C)%*%Compute.Q(Cov,return.proj=T)%*%out$C)
  ind.tmp <- which(!is.na(out$L[,1]))
  eigs <- eigen(1/p*R%*%t(out$L[ind.tmp,])%*%out$L[ind.tmp,]%*%t(R),symmetric=T)
  out$stand.eigen.values <- (n-d)*eigs$values/mean(out$Sigma,na.rm=T)
  out$Perc.Var.Explained <- out$eigen.values/(sum(out$eigen.values) + mean(out$Sigma,na.rm=T))
  M <- t(R)%*%eigs$vectors
  out$L <- out$L%*%M
  out$C <- out$C%*%solve(t(M))
  out$Var.L[ind.tmp] <- lapply(out$Var.L[ind.tmp],function(x){t(M)%*%x%*%M})
  return(out)
}