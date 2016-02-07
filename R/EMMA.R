EMMA <- function (y, X = NULL, Z = NULL, K = NULL, REML = TRUE, silent = FALSE){
  if (!silent) {
    count <- 0
    tot <- 1
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  y.or <- y
  if (is.null(X)) {
    x.or <- matrix(rep(1, length(y)), ncol = 1)
  }
  else {
    x.or <- X
  }
  if (is.null(Z)) {
    z.or <- diag(length(y))
  }
  else {
    z.or <- Z
  }
  good <- which(!is.na(y))
  y <- y[good]
  Z <- Z[good, ]
  if (is.null(X)) {
    X <- matrix(rep(1, length(y)), ncol = 1)
  }
  else {
    X <- as.matrix(X[good, ])
  }
  if (!is.null(X)) {
    if (is.list(X)) {
      if (is.list(X[[1]])) {
        X = X[[1]][[1]]
      }
      else {
        X = X[[1]]
      }
    }
    else {
      X = X
    }
  }
  if (is.null(Z)) {
    Z <- diag(length(y))
  }
  if (!is.null(Z)) {
    if (is.list(Z)) {
      if (is.list(Z[[1]])) {
        Z = Z[[1]][[1]]
      }
      else {
        Z = Z[[1]]
      }
    }
    else {
      Z = Z
    }
  }
  if (is.null(K)) {
    K <- diag(length(y))
  }
  if (!is.null(K)) {
    if (is.list(K)) {
      if (is.list(K[[1]])) {
        K = K[[1]][[1]]
      }
      else {
        K = K[[1]]
      }
    }
    else {
      K = K
    }
  }
  n = length(y)
  p <- dim(X)[2]
  II <- as(diag(n), "sparseMatrix")
  S <- diag(n) - X %*% crossprod(solve(crossprod(X)), t(X))
  S <- as(S, "sparseMatrix")
  K <- as(K, "sparseMatrix")
  Z <- as(Z, "sparseMatrix")
  de <- log(n)
  H <- (Z %*% (K %*% t(Z))) + (de * II)
  eigSHS <- eigen(S %*% H %*% S, symmetric = TRUE)
  lambda <- eigSHS$values[1:(n - p)] - de
  U <- eigSHS$vectors[, (1:(n - p))]
  Hb.system <- eigen(H, symmetric = TRUE)
  phi <- Hb.system$values - de
  eta <- crossprod(U, y)
  if (REML == TRUE) {
    minimfunc <- function(delta) {
      (n - p) * log(sum(eta^2/(lambda + delta))) + sum(log(lambda + 
                                                             delta))
    }
  }
  else {
    minimfunc <- function(delta) {
      n * log(sum(eta^2/(lambda + delta))) + sum(log(phi + 
                                                       delta))
    }
  }
  REML <- optimize(minimfunc, lower = 9^(-9), upper = 9^9, 
                   tol = 1e-06)
  delta <- REML$minimum
  hhhh <- Z %*% K %*% t(Z) + delta * diag(length(y))
  H.hat.inv <- solve(as(hhhh, "sparseMatrix"))
  aaa <- t(X) %*% H.hat.inv %*% X
  bbb <- t(X) %*% H.hat.inv %*% y
  beta <- solve(aaa, bbb)
  error <- y - (X %*% beta)
  sigma2.u <- sum((eta^2)/(lambda + delta))/(n - p)
  sigma2.e <- delta * sigma2.u
  u <- t(Z %*% K) %*% (H.hat.inv %*% error)
  df <- n - p
  pi <- 3.14159
  ll <- -0.5 * (REML$objective + df + (df * log((2 * pi)/df)))
  Vi <- (1/sigma2.u) * H.hat.inv
  P <- Vi - ((Vi %*% X) %*% solve(crossprod(X, Vi %*% X), crossprod(X, 
                                                                    Vi)))
  Var.u <- (sigma2.u^2) * (crossprod(Z %*% K, P) %*% (Z %*% 
                                                        K))
  PEV.u <- sigma2.u * K - Var.u
  var.beta <- solve(crossprod(X, Vi %*% X))
  out <- matrix(c(sigma2.u, sigma2.e))
  rownames(out) <- c("V(u)", "V(e)")
  AIC = (-2 * ll) + (2 * dim(X)[2])
  BIC = (-2 * ll) + (log(length(y)) * dim(X)[2])
  fitted.y <- x.or %*% beta
  fitted.y <- fitted.y + (z.or %*% u)
  fitted.u <- (z.or %*% u)
  fitted.y.good <- fitted.y[good]
  residuals2 <- y.or[good] - fitted.y[good]
  rownames(u) <- colnames(Z)
  rownames(beta) <- colnames(X)
  if (!silent) {
    count <- count + 1
    setTxtProgressBar(pb, (tot/tot))
  }
  return(list(var.comp = out, beta.hat = beta, u.hat = u, Var.u.hat = (Var.u), 
              Var.beta.hat = var.beta, PEV.u.hat = (PEV.u), LL = ll, 
              AIC = AIC, BIC = BIC, V.inv = H.hat.inv, X = X, Z = Z, 
              K = K, fitted.y = fitted.y, fitted.u = fitted.u, residuals = error, 
              cond.residuals = residuals2, fitted.y.good = fitted.y.good))
}