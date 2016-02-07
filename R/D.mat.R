D.mat <- function (X, min.MAF = NULL, max.missing = NULL, impute.method = "mean", 
          tol = 0.02, n.core = 1, shrink = FALSE, return.imputed = FALSE, 
          ploidy = 2){
  ty <- apply(X, 2, function(x) {
    length(table(x))
  })
  vv <- which(ty > 2)
  X2 <- X[, vv]
  if (ploidy == 2) {
    X3 <- apply(X2, 2, function(x) {
      y <- x
      y[which(x == 1 | x == -1)] <- 0
      y[which(x == 0)] <- 1
      return(y)
    })
    X6 <- A.mat(X3, min.MAF = min.MAF, max.missing = max.missing, 
                impute.method = impute.method, tol = tol, n.core = n.core, 
                shrink = shrink, return.imputed = return.imputed)
  }
  else {
    X3 <- X2 - (ploidy/2)
    possible <- (-(ploidy/2):(ploidy/2))
    homo <- c(possible[1], possible[length(possible)])
    hete <- setdiff(possible, homo)
    X4 <- apply(X3, 2, function(x) {
      y <- x
      y[which(x %in% hete)] <- 1
      y[which(x %in% homo)] <- 0
      return(y)
    })
    ty2 <- apply(X4, 2, function(x) {
      length(table(x))
    })
    vv2 <- which(ty2 > 1)
    X5 <- X4[, vv2]
    X6 <- A.mat(X5, min.MAF = min.MAF, max.missing = max.missing, 
                impute.method = impute.method, tol = tol, n.core = n.core, 
                shrink = shrink, return.imputed = return.imputed)
  }
  return(X6)
}