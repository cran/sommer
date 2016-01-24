

score.calc <- function(M, y, ZO, Hinv, min.MAF, X2, P3D=TRUE, method="AI", iters=50, R=NULL, REML=TRUE, draw=TRUE) {
  #y <- y[which(!is.na(y))]
  #M <- M[which(!is.na(y)),]
  #ZO <- ZO[which(!is.na(y)),]
  #X2 <- X2[which(!is.na(y)),]
  Z <- diag(dim(M)[1])
  # Z is square diagonal for n inds, M are the marker matrix
  dim(M); dim(Z)
  m <- ncol(M); m # number of markers
  scores <- array(0, m) # number of scores required for each marker
  ##########################
  ##########################
  # for each marker do this
  ##########################
  ##########################
  ####################
  ## initialize the progress bar
  count <- 0
  tot <- m
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  #####################
  for (i in 1:m) {
    ###################
    count <- count + 1
    ###################
    Mi <- M[, i] ; Mi# extract the vector or column of the ith marker
    freq <- mean(Mi + 1, na.rm = TRUE)/2
    MAF <- min(freq, 1 - freq); MAF #MAF of such marker
    if (MAF >= min.MAF) {
      not.NA.gid <- which(!is.na(Mi)) # with no missing data
      temp <- rep(1, length(Mi)) # temporal vector, 0 if has cas data, 1 if is NA 
      temp[not.NA.gid] <- 0
      not.NA.obs <- which(Z %*% temp != 1) # if are != 1 has actual data 
      n2 <- length(not.NA.obs); n2 # #inds with data 
      y2 <- matrix(y[not.NA.obs], n2, 1); y2 # y in matrix form
      Z2 <- Z[not.NA.obs, not.NA.gid] ; dim(Z2)# Z with data
      # bind X (intercept)  with Z Mi product of Z and calls in Mi marker
      # [X ZMi] where Mi is the vector for each marker analyzed
      X3 <- cbind(X2[not.NA.obs, ], Z2 %*% Mi[not.NA.gid]); dim(X3); X3[1:10,]
      
      p <- ncol(X3); p # number of columns
      v1 <- 1
      v2 <- n2 - p; v2 # number of inds - numb of columns = degrees of freedom
      ############################
      
      if (P3D) { #if P3D, meaning normal model, variance components estimated once
        H2inv <- Hinv[not.NA.obs, not.NA.obs]
      }else{ # if user wants to estimate variance components for each marker
        
        if(method == "EMMMA"){
          random <- which(unlist(lapply(ZO, function(x){names(x)[1]})) == "Z") # elements of Z that are RANDOM
          if(length(random) > 0){# EMMA buth in 2-level list provided
            H2inv  <- EMMA(y=y, X=X3, Z=ZO[[random]][[1]], K=ZO[[random]][[2]],  REML=REML)$V.inv 
          }else{
            H2inv <- EMMA(y=y, X=X3, Z=ZO[[1]], K=ZO[[2]],REML=REML)$V.inv 
          }
        }
        if(method == "EM"){
          H2inv <- EM(y=y, X=X3, ETA=ZO, iters = iters, REML=REML, draw=draw, silent=TRUE)$V.inv
        }
        if(method == "AI"){
          H2inv <- AI(y=y, X=X3, ZETA=ZO, R=R, REML=REML, draw=draw, silent=TRUE)$V.inv
        }
      }
      
      ############################
      # XH-X'
      W <- crossprod(X3, H2inv %*% X3)
      # [XH-X']-
      Winv <- try(solve(W), silent = TRUE)
      if (class(Winv) != "try-error") {
        # [XH-X']- XH-y = Beta
        beta <- Winv %*% crossprod(X3, H2inv %*% y2)
        # e = y - XB
        resid <- y2 - X3 %*% beta
        # e H- e / n-p = SSe/(n-p)
        s2 <- as.double(crossprod(resid, H2inv %*% 
                                    resid))/v2
        # Variance covariance SSe/(n-p) * [XH-X']-
        CovBeta <- s2 * Winv
        # F statistic is Beta^2 / Var(Beta)
        Fstat <- beta[p]^2/CovBeta[p, p]
        # (n-p) / (n-p + 1 * F)
        x <- v2/(v2 + v1 * Fstat)
        # -log10 of beta distribution
        scores[i] <- -log10(pbeta(x, v2/2, v1/2))
      }
    }
    ################################
    setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    ################################
  }
  return(scores)
}