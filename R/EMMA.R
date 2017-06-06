EMMA <- function (y, X=NULL, ZETA=NULL, REML=TRUE, silent=FALSE, che=TRUE, forced=NULL, EIGEND=FALSE){
  
  NNN <- length(ZETA)
  if(EIGEND){
    if(NNN == 1){
      DISO <- dim(ZETA[[1]]$Z)
      if(DISO[1] != DISO[2]){
        stop("EIGEN DECOMPOSITION EIGEND ONLY WORKS FOR SQUARE PROBLEMS 
             'Z' MATRIX IS NOT SQUARE",call.=FALSE)
        
      } 
    }else{
      stop("EIGEND FEATURE ONLY WORKS FOR ONE RANDOM EFFECT",call.=FALSE)
      jkl <- c(23,18,9,20,20,5,14, NA,2,25,NA,7,9,15,22,1,14,14,25,NA,3,15,22,1,18,18,21,2,9,1,19)
      oh.yeah <- paste(letters[jkl],collapse = "")
    }
  }
  
  y.or <- y
  
  ### make full function
  '%!in%' <- function(x,y)!('%in%'(x,y))
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  ## MTG chose
  if(EIGEND){
    cat("EIGEND feature activated. Eigen decomposition of K will be performed\n")
  }
  ###########################
  ## y is a vector for the response variable
  ## X is an incidence matrix for fixed effects
  ## Z is a list of lists for each random effect
  # the list of Z can or cannot include the covariance matrix for such random effect
  # if provided must be provided as Z=list(list(Z=Z1,K=K1),list(Z=Z2,K=K2), etc) 
  ############################
#   if(che){ # if coming from mmer don't check
#     if(is.list(ZETA)){
#       if(is.list(ZETA[[1]])){ # if was provided as a two level list
#         ZETA=ZETA
#       }else{ # if was provided as a one level list
#         ZETA=list(ZETA)
#       }
#     }else{
#       #stop;
#       cat("\nThe random effects need to be provided in a list format, please see examples")
#     }
#   }
  ###########################
  # if X matrix is not present
  if(is.null(X) & is.null(ZETA)){ # nothing in the model
    tn = length(y); xm <- matrix(1,tn,1)
    yv <- scale(y)
    res <- lm(yv~xm-1) # intercept model
  }else{ ### GOOD AND NORMAL MODEL ###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(is.null(X) & !is.null(ZETA)){ # only random effects present
      tn = length(y); X <- matrix(1,tn,1)
    }
    if(!is.null(X) & !is.null(ZETA)){ # both present, extract xm from X, check double list
      if(is.list(X)){
        if(is.list(X[[1]])){
          X=X[[1]][[1]]
        }else{
          X=X[[1]]
        }
      }else{
        X=as.matrix(X) 
      }
    }
    ############################################
    ## if K matrices are not present in ZETA
    # add an identity matrix to all effects in ETA that did not provide a var-cov matrix
    #if(is.null(R)){R <- diag(length(y))} # dim(x[[1]])[2]
    
    if(che){ # if needs to be checked fill empty matrices, else just skip
      ZETA <- lapply(ZETA, function(x){
        if(length(x) == 1){
          provided <- names(x)
          if(provided == "Z"){
            y <- list(Z=x[[1]],K=diag(dim(x[[1]])[2]))
          }
          if(provided == "K"){
            y <- list(Z=diag(length(y)),K=x[[1]])
          }else{
            stop("Names of matrices provided can only be 'Z' or 'K', the names you provided don't match the arguments required",call.=FALSE)
            
          }
        }else{y <- x}; 
        return(y)
      })
    }
    #######################################################
    ## order random effects according to degrees of freedom
    tokeep <- names(ZETA)
    df <- unlist(lapply(ZETA, function(x){dim(x[[1]])[2]}))
    df2 <- sort(df, decreasing = FALSE)
    df.ord <- numeric() # # 
    for(u in 1:length(df)){
      df.ord[u] <- which(df2 %in% df[u])[1]
      df2[df.ord[u]] <- NA
    }
    ZETA <- ZETA[df.ord]
    names(ZETA) <- tokeep[df.ord]
    #nano <- tokeep[df.ord]
    #print(tokeep[df.ord])
    if(!is.null(forced)){forced <- forced[c(df.ord, (length(forced)))]}
    #####################################################
    ## to use later for fitted values
    X.or <- as.matrix(X)
    zeta.or <- ZETA
    zeta.or  <- lapply(zeta.or , function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    ## if we have eigen structure impute response
    if((length(ZETA)==1) & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2]) & EIGEND){
      misso <- which(is.na(y))
      if(length(misso) >0){
        y[misso] <- median(y, na.rm=TRUE)
      }
    }
    ZETA2 <- ZETA; y2 <- y ; good <- which(!is.na(y)) # make a copy and identify no missing data
    #ZETA <- lapply(ZETA2, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
    if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
      
      ZETA <- lapply(ZETA2, function(x,good){
        if(dim(x[[1]])[2] == dim(x[[2]])[2]){ # if square
          x[[1]] <- x[[1]][good,good]; x[[2]]<- x[[2]][good,good]
          return(x)
        }else{ # if general mixed model
          x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]
          return(x)
        }}, good=good)
    }else{
      ZETA <- lapply(ZETA2, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
    }
    ################
    y <- y[good]
    ZETA <- lapply(ZETA, function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    X <- as.matrix(X[good,])
    Xt <- t(X)
    #R <- R[good,good]
    NNN <- length(ZETA)
    
    ################################################################################################
    ################################################################################################
    ############ START ALGORITHM
    ################################################################################################
    ################################################################################################
    if(!silent){
      count <- 0
      tot <- 3
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    
    
    NNN <- length(ZETA)
    Zlist <- lapply(ZETA, function(x){x$Z})
    Klist <- lapply(ZETA, function(x){x$K})
    
    if(EIGEND){
      if(NNN == 1){
        # eigen decomposition doesn't take missing values so we impute
        missed <- which(is.na(y))
        if(length(missed)>0){y[missed] <- mean(y, na.rm=TRUE)}
        K <- do.call("adiag1",Klist)
        ED <- eigen(K)
        Us <- ED$vectors
        Ds <- diag(ED$values)
        
        yyy <- y
        xxx <- X
        kkk <- K
        y <- as.vector(t(Us) %*% as.matrix(y,ncol=1))
        X <- t(Us) %*% X
        K <- Ds 
      }else{
        stop("EIGEND FEATURE ONLY WORKS FOR ONE RANDOM EFFECT",call.=FALSE)
      }
    }
    
    q = dim(X)[2]
    n = length(y)
    lz <- length(Zlist)
    spI <- diag(n)
    S <- spI - tcrossprod(X %*% solve(crossprod(X)), X) ## S = I - X(X'X)-X'
    #######################################################################
    ############  START MULTIPLE RANDOM EFFECTS AND FORCED ################
    if(is.null(forced) & NNN > 1){  
      Z <- do.call("cbind", Zlist)
      ####========MIN.FUNCTION===========####
      minimfunctionouter <- function(weights = rep(1/lz, lz)) {
        weights = weights/sum(weights)
        ZKZt <- matrix(0, nrow = n, ncol = n)
        for (i in 1:lz) {
          ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]] %*% 
                                                   Klist[[i]], Zlist[[i]])
        }
        offset <- log(n)
        ZKZtandoffset <- ZKZt + offset * spI
        SZKZtSandoffset <- {
          S %*% ZKZtandoffset
        } %*% S
        svdSZKZtSandspI <- eigen(SZKZtSandoffset, symmetric = TRUE)
        Ur <- svdSZKZtSandspI$vectors[, 1:(n - q)]
        lambda <- svdSZKZtSandspI$values[1:(n - q)] - offset
        eta <- crossprod(Ur, y)
        minimfunc <- function(delta) {
          (n - q) * log(sum(eta^2/{
            lambda + delta
          })) + sum(log(lambda + delta))
        }
        optimout <- optimize(minimfunc, lower = 0, upper = 10000, 
                             tol = 1e-09)
        return(optimout$objective)
      }####========MIN.FUNCTION===========####
      
      
      weights <- optim(par = rep(1/lz, lz), fn = minimfunctionouter, 
                       method = "L-BFGS-B", lower = rep(0, lz), upper = rep(1,lz))$par
      
      #weights <- weights#/length(weights)
      
      zero.comp <- which(weights==0)
      no.zero.comp <- which(weights!=0)
      lz2 <- length(no.zero.comp)
      if(length(zero.comp)>0){
        cat("\nOne or more variance components pushing to be zero. Boundary constraint applied.\n")
        ####========MIN.FUNCTION===========####
        minimfunctionouter2 <- function(weights = rep(1/lz2, lz2)) {
          weights = weights/sum(weights)
          ZKZt <- matrix(0, nrow = n, ncol = n)
          for (i in no.zero.comp) {
            ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]] %*% Klist[[i]], Zlist[[i]])
          }
          offset <- log(n)
          ZKZtandoffset <- ZKZt + offset * spI
          SZKZtSandoffset <- {
            S %*% ZKZtandoffset
          } %*% S
          svdSZKZtSandspI <- eigen(SZKZtSandoffset, symmetric = TRUE)
          Ur <- svdSZKZtSandspI$vectors[, 1:(n - q)]
          lambda <- svdSZKZtSandspI$values[1:(n - q)] - offset
          eta <- crossprod(Ur, y)
          minimfunc <- function(delta) {
            (n - q) * log(sum(eta^2/{lambda + delta})) + sum(log(lambda + delta))
          }
          optimout <- optimize(minimfunc, lower = 0, upper = 10000, 
                               tol = 1e-09)
          return(optimout$objective)
        }####========MIN.FUNCTION===========####
        weights2 <- optim(par = rep(1/lz2, lz2), fn = minimfunctionouter2, 
                          method = "L-BFGS-B", lower = rep(0, lz2), upper = rep(1,lz2))$par
        weights[zero.comp] <- 0
        
        ## when multiple random effects and some var.comp are zero
        #if((length(weights2)==1) & (weights2 == 1)){
        #  print("oh")
        #  dddd <- EMMA(y, X=X, ZETA=ZETA[-zero.comp])$var.comp
        #  weights[no.zero.comp] <- dddd[-length(dddd),]
        #}else{
          weights[no.zero.comp] <- weights2#/length(weights)
        #}
        
      }# end of variance components being zero
      
      #if(!silent){
      #  setTxtProgressBar(pb, (1/tot))### keep filling the progress bar
      #}
      
      sigma <- var(y) * weights
    }else{ ## IF FORCED
      weights <- 1
      if(!is.null(forced)){
        weights <- forced[1:length(ZETA)]
        sigma <- var(y) * weights
      }
      
      #if(!silent){
      #  setTxtProgressBar(pb, (2/tot))### keep filling the progress bar
      #}
      
    }### END OF FORCING
    ##############  END MULTIPLE RANDOM EFFECTS AND FORCED ################
    #######################################################################
    
    ##### once var.comp have been obtained create V matrix
    Z <- do.call("cbind", Zlist)
    
    weights <- weights/sum(weights)
    ZKZt <- matrix(0, nrow = n, ncol = n)
    Klistweighted <- Klist
    ### weights are var.comp
    ### form ZKZ by summing each Zi.Ki.Zi
    for (i in 1:lz) {
      Klistweighted[[i]] <- weights[i] * Klist[[i]]
      ZKZt <- ZKZt + weights[i] * tcrossprod(Zlist[[i]] %*% 
                                               Klist[[i]], Zlist[[i]])
    }
    ## bind in a diagonal all K's*sigma
    K <- do.call("adiag1",Klistweighted) # .bdiag(Klistweighted)
    ZK <- as.matrix(Z %*% K) ## ZK
    offset <- log(n)
    
    ## -----------------------------------------------------------------##
    ##-------------------- V = ZKZ' + delta.prov*I ---------------------##
    ## -----------------------------------------------------------------##
    ZKZtandoffset <- ZKZt + offset * spI ## ZKZ' + I*offset
    
    ## S = I - X(X'X)-X' and then
    ## S ZKZ' + I*offset S
    SZKZtSandoffset <- {S %*% ZKZtandoffset } %*% S
    ## -----------------------------------------------------------##
    ## ------------------- eigen(SHS) --> UV ---------------------##
    ## -----------------------------------------------------------##
    svdSZKZtSandspI <- eigen(SZKZtSandoffset, symmetric = TRUE)
    
    Ur <- svdSZKZtSandspI$vectors[, 1:(n - q)] # eigen vectors of SHS = [I - X(X'X)-X'] [ZKZ' + I*theta] [I - X(X'X)-X']
    lambda <- svdSZKZtSandspI$values[1:(n - q)] - offset # eigen values of SHS = [I - X(X'X)-X'] [ZKZ' + I*theta] [I - X(X'X)-X']
    
    # if ML
    if(!REML){
      Hb.system <- eigen(SZKZtSandoffset, symmetric = TRUE)
      phi <- Hb.system$values - offset
    }
    
    ## -----------------------------------------------------------##
    ##-------------------- estimate parameter eta and delta ---------------------## 
    ## part of "y" variability attributed to random effects 
    ## -----------------------------------------------------------##
    eta <- crossprod(Ur, y)
    
    ## -----------------------------------------------------------##
    ##-------------------- likelihood function using eta, lambda and delta ---------------------##
    ## -----------------------------------------------------------##
    if(REML){
      minimfunc <- function(delta) {(n - q) * log(sum(eta^2/{lambda + delta})) + sum(log(lambda + delta)) }
    }else{
      minimfunc <- function(delta) {n * log(sum(eta^2/(lambda + delta))) + sum(log(phi + delta)) }
    }
    
    ## ------------------------------------------------------------------------------##
    ##-------------------- REML estimator of delta=Var(e)/Var(u)---------------------##
    ## ------------------------------------------------------------------------------##
    optimout <- optimize(minimfunc, lower = 0, upper = 10000, 
                         tol = 1e-09)
    deltahat <- optimout$minimum
    
    #if(!silent){
    #  setTxtProgressBar(pb, (2/tot))### keep filling the progress bar
    #}
    ##-------------------- H- = (H = ZKZ' + delta*I)- ---------------------## 
    # remember V- = 1/var(g) * H-  which solves the random effect model, here there's an extra ridge parameter
    Hinvhat <- solve(ZKZt + deltahat * spI)
    XtHinvhat <- crossprod(X, Hinvhat)
    betahat <- solve(XtHinvhat %*% X, XtHinvhat %*% y)
    ## var beta
    var.beta <- solve(crossprod(X, Hinvhat %*% X))
    ## residuals
    ehat <- (y - {X %*% betahat})
    ## V-(y-Xb)
    Hinvhatehat <- Hinvhat %*% ehat
    ## sigma u
    sigmausqhat <- sum(eta^2/{lambda + deltahat})/(n - q)
    ## sigma error
    sigmaesqhat <- deltahat * sigmausqhat
    
    if(NNN > 1){
      sigma <- c(sigma,sigmaesqhat)
    }else{
      sigma <- c(sigmausqhat,sigmaesqhat)
    }
    
    ## blups u
    if(EIGEND){
      Usi <- solve(t(Us))
      uhat <- Usi %*% crossprod(ZK, Hinvhatehat)
    }else{
      #u <- t( Z %*% K ) %*% (H.hat.inv %*% error)  # ZK H- (y-XB) # gianola uses letter V- instead of H-
      uhat <- crossprod(ZK, Hinvhatehat)
    }
    
    
    
    Vinv <- (1/sigmausqhat) * Hinvhat
    namesuhat <- c()
    for (i in 1:length(Klist)) {
      namesuhat <- c(namesuhat, paste(paste("K", i, sep = "_"),colnames(Klist[[i]]), sep = ""))
    }
    ## DEGREES OF FREEDOM
    df <- n - q
    ## LOG-LIKELIHOOD
    loglik <- -0.5 * (optimout$objective + df + df * log(2 * pi/df))
    ## AIC and BIC
    AIC = (-2 * loglik ) + ( 2 * dim(X)[2]) # k=2, npar=2
    BIC = (-2 * loglik ) + ( log(length(y)) * dim(X)[2])
    ## VAR(U.hat) and PEV(u.hat)
    
    if(EIGEND){ # sigma^4  Ui' (ZKP ZK) Ui
      P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
      Var.u <-   (sigmausqhat^2) *  Usi  %*% tcrossprod( crossprod(ZK, P)  %*%  (ZK),Usi  )
      PEV.u <-  sigmausqhat * K - varuhat  # standard errors (SE) for each individual
    }else{
      P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv))
      varuhat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK  # sigma^4 ZKP ZK
      PEVuhat <- sigmausqhat * K - varuhat  # standard errors (SE) for each individual
    }
    
    ### Var(b.hat)
    varbetahat <- solve(crossprod(X, Vinv %*% X))
    ### EXTRAS
    fitted.y <- X.or %*% betahat
    zz.list <- lapply(zeta.or,function(x){x$Z})
    z.or <- do.call("cbind",zz.list)
    fitted.y <- fitted.y + (z.or %*% uhat)
    fitted.u <- (z.or %*% uhat)
    fitted.y.good <- fitted.y[good]
    residuals2 <- y.or[good] - fitted.y[good] # conditional residuals
    rownames(uhat) <- colnames(Z)
    rownames(betahat) <- colnames(X)
    
    if(!silent){
      setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
    }
    
    #### tests
    #if (test) {
    #  Xsqtestu <- uhat^2/diag(varuhat)
    #  puhat <- pchisq(Xsqtestu, df = 1, lower.tail = F, log.p = F)
    #  p.adjust.M <- p.adjust.methods
    #  p.adjuhat <- sapply(p.adjust.M, function(meth) p.adjust(puhat, 
    #                                                          meth))
    #  Xsqtestbeta <- betahat^2/diag(varbetahat)
    #  pbetahat <- pchisq(Xsqtestbeta, df = 1, lower.tail = F, 
    #                     log.p = F)
    #  p.adjbetahat <- sapply(p.adjust.M, function(meth) p.adjust(pbetahat, meth))
    #}
    uhat <- as.numeric(uhat)
    names(uhat) <- namesuhat
    
    u.hatl <- list()
    varuhatl <- list()
    PEVuhatl <- list()
    levo <- c(unlist(lapply(ZETA,function(x){dim(x$Z)[2]})))
    for(ss in 1:length(levo)){
      u.hatl[[ss]] <- as.matrix(uhat[1:levo[ss]])
      rownames(u.hatl[[ss]]) <- colnames(ZETA[[ss]]$Z)
      ## var.uhat
      varuhatl[[ss]] <- varuhat[1:levo[ss],1:levo[ss]]
      rownames(varuhatl[[ss]]) <- colnames(ZETA[[ss]]$Z)
      colnames(varuhatl[[ss]]) <- rownames(varuhatl[[ss]])
      ##pev
      PEVuhatl[[ss]] <- PEVuhat[1:levo[ss],1:levo[ss]]
      rownames(PEVuhatl[[ss]]) <- colnames(ZETA[[ss]]$Z)
      colnames(PEVuhatl[[ss]]) <- rownames(PEVuhatl[[ss]])
      ## clean data
      varuhat <- varuhat[-c(1:levo[ss]),-c(1:levo[ss])]
      PEVuhat <- PEVuhat[-c(1:levo[ss]),-c(1:levo[ss])]
      uhat <- uhat[-c(1:levo[ss])]
    }
    #print(tokeep[df.ord])
    #print(length(u.hatl))
    names(u.hatl) <- c(tokeep[df.ord])
    #print(names(u.hatl))
    
    if(EIGEND){
      ehat <- Usi %*% ehat
      #H.hat.inv <- Usi  %*% H.hat.inv %*% t(Usi)
      H.hat.inv <- Usi%*%tcrossprod(Hinvhat, Usi)
    }
    sigma <- matrix(sigma)
    rownames(sigma) <- c(paste("V(u.",1:length(ZETA),")",sep=""),"V(Error)")
    
    out1 <- as.data.frame(sigma)
    out1[,"constraint"] <- "Positive"
    convergence<-TRUE
    res <- list(var.comp=out1, beta.hat = betahat, 
                u.hat = u.hatl, Var.u.hat = varuhatl, 
                Var.beta.hat = var.beta, PEV.u.hat = PEVuhatl, 
                LL = loglik, AIC=AIC, BIC=BIC, V.inv=Hinvhat, X=X, Z=Z, K=K, 
                fitted.y=fitted.y, fitted.u=fitted.u, residuals=ehat, 
                cond.residuals=residuals2, fitted.y.good=fitted.y.good,
                convergence=convergence)
  } ### end of model with random effects
  return(res)
}