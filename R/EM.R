EM <- function(y,X=NULL,ZETA=NULL,R=NULL,iters=30,draw=TRUE,silent=FALSE, constraint=TRUE, init=NULL, forced=NULL, tolpar = 1e-04, tolparinv = 1e-06){
  convergence <- FALSE
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### fix that some initially were not specified
  if(is.null(X)){
    X <- matrix(1,nrow=length(y))
  }
  if(is.null(R)){
    R <- list(units=diag(length(y)))
  }
  y.or <- y
  X.or <- X
  ZETA.or <- ZETA
  R.or <- R
  if(is.null(names(ZETA))){
    varosss <- c(paste("u.",1:length(ZETA), sep=""))
  }else{
    varosss <- c(names(ZETA))
  }; varosssZ <- varosss
  if(is.null(names(R))){
    varosss <- c(varosss,paste("Res.",1:length(R),sep=""))
  }else{
    varosss <- c(varosss,names(R))
  }
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  
  #### now reduce the original inputs
  #y <- scale(y)
  good <- which(!is.na(y))
  y <- y[good]
  X <- as.matrix(X[good,])
  ZETA <- lapply(ZETA, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
  R <- lapply(R, function(x,good){x <- x[good,good]; return(x)}, good=good)
  #### get sizes of reduced
  qr <- qr(X)
  ranx <- length(y)-qr$rank # length(good)
  nz <- length(ZETA)
  nr <- length(R)
  nx <- dim(X)[2]
  dimzs <- lapply(ZETA,function(x){dim(as.matrix(x$Z))[2]})
  dimrs <- lapply(R,function(x){dim(x)[1]})
  N <- length(y)
  
  #### get the indexes for all effects
  sta <- 1
  ind1 <- numeric()
  for(u in 1:length(dimzs)){
    sta <- dimzs[[u]]+sta
    ind1[u] <- sta
  }; ind1<- c(1,ind1[-c(length(ind1))]);
  ind2 <- (c(ind1-1, sum(unlist(dimzs))))[-1]
  
  ind1g <- ind1 + nx ## indexes for Gi
  ind2g <- ind2 + nx ## indexes for Gi
  ind1 <- c(1,ind1 + nx) ## indexes for all including x
  ind2 <- c(nx,ind2 + nx) ## indexes for all including x
  
  #### initialize var components
  if(is.null(init)){
    varcom <- rep(var(y,na.rm=TRUE)/(nz+nr),nz+nr)
    names(varcom) <- c(rep("z",nz),rep("r",nr))
    varcomz <- varcom[which(names(varcom)=="z")]
    varcomr <- varcom[which(names(varcom)=="r")]
  }else{
    varcom <- init 
    if(length(init)!= (nr+nz)){
      stop("Please provide initial values for all variance components",call. = FALSE)
    }else{
      names(varcom) <- c(rep("z",nz),rep("r",nr))
      varcomz <- varcom[which(names(varcom)=="z")]
      varcomr <- varcom[which(names(varcom)=="r")]
    }
  }
  
  ll2=-10000000 # initial log-likelihood
  ll.stored <- ll2
  conv=0 # convergence
  wi=0 # counter
  taper <- rep(1, iters) # weighting parameter for updates
  #taper[1:3] <- c(.5,.7,1)#c(0.5, 0.7) # c(0.5, 0.7)
  
  #### get inverse of covariance and residual matrices
  Riw <-  lapply(R,function(x){solve(x)})
  Ki <- lapply(ZETA,function(x){
    findrank <- qr(x$K)$rank
    if(findrank < dim(x$K)[1]){
      return(solve(x$K + diag(tolparinv,dim(x$K)[1])))
    }else{
      return(solve(x$K))
    }
  })
  
  varso <- as.matrix(varcom)
  #### start
  if(!silent){
    count <- 0
    tot <- iters+1
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  if(!is.null(forced)){
    varcom <- forced
    if(length(init)!= (nr+nz)){
      stop("Please provide initial values for all variance components",call. = FALSE)
    }else{
      names(varcom) <- c(rep("z",nz),rep("r",nr))
      varcomz <- varcom[which(names(varcom)=="z")]
      varcomr <- varcom[which(names(varcom)=="r")]
    }
    iters=1
  }
  ################### =================================================================
  ################### =================================================================
  ################### =================================================================
  while (conv==0) { # ==================== START EM ALGORITHM =========================
    wi=wi+1
    if(!silent){
      count <- count + 1
    }
    ##############################
    ### for R * se as a direct sum
    Rse <- R[[1]]*0
    for(u in 1:length(R)){
      Rse <- Rse + R[[u]]
    }
    ##############################
    ## R inverse
    Rsei <- solve(Rse)#solve(Rse)
    ##############################
    ## G inverse
    varoz <- as.list(varcomz)  # variance components as list, no error included
    Gi <- lapply(as.list(c(1:nz)),function(x,Ko,v){
      oo=Ko[[x]]*as.numeric(varcomr[1]/(v[[x]])); return(oo)
    }, Ko=Ki, v=varoz) ## K*v(u)
    
    Zbind <- do.call(cbind,lapply(ZETA, function(x){x$Z}))
    XZbind <- as(cbind(X,Zbind), Class="sparseMatrix")
    CM <- t(XZbind) %*% Rsei %*% XZbind
    
    ## add Gi's to the C matrix
    for(u in 1:length(Gi)){
      rox <- ind1g[u]; cox <- ind2g[u]
      CM[rox:cox,rox:cox] <- CM[rox:cox,rox:cox] + Gi[[u]]
    }
    
    ## do gaussian eliminations (absorption)
    RHS <- t(XZbind) %*% Rsei %*% y
    ge <- solve(CM,RHS)
    ## invert the coefficient matrix
    CMi <- solve(CM)
    
    ##%%%%%%%%%%%%%%
    ### likelihood
    ##%%%%%%%%%%%%%%
    na <- sum(unlist(dimzs));na
    se <- varcomr
    #### Calculate ypy by building M
    left <- t(XZbind) %*% Rsei %*% y
    top <- t(left)
    corner <- t(y)%*%Rsei%*%y
    
    M <- rbind(cbind(corner,top),cbind(left,CM))
    M[1:4,1:4]
    vb <- try(chol(as(M, Class="sparseMatrix"),pivot=FALSE), silent=TRUE)
    if(class(vb)=="try-error"){
      vb <- try(chol(as(M+diag(tolparinv,dim(M)[1]), Class="sparseMatrix")), silent=TRUE)
      if(class(vb)=="try-error"){
        ypy <- 1#(L[1,1]^2); ypy
        logdc <- 2#2* sum(log(diag(L)[-1])); logdc
      }else{
        L <- t(vb); L[1:4,1:4]
        ypy <- (L[1,1]^2); ypy
        logdc <- 2* sum(log(diag(L)[-1])); logdc
      }
    }else{
      L <- t(vb); L[1:4,1:4]
      ypy <- (L[1,1]^2); ypy
      logdc <- 2* sum(log(diag(L)[-1])); logdc
    }
    
    
    logda <- sum(unlist(lapply(ZETA,function(x){determinant(x$K, logarithm = TRUE)$modulus[[1]]})))
    nlogsu <- log(varcomz)*unlist(dimzs) # n*log(su)
    
    
    ll <- - 0.5 * ( ((N - ranx - na)*log(se)) + logdc + logda + sum(nlogsu) +  ypy);ll
    if (abs(ll-ll2) < tolpar | wi == iters ){ ## CONVERGENCE, not sure if should be absolute value or not
      conv=1
      if (abs(ll-ll2) < tolpar){
        convergence <- TRUE
      }
    }
    ll2=(ll) # replace the initial logLik value for the value reached
    ll.stored <- c(ll.stored,ll2)
    ##%%%%%%%%%%%%%%
    ##%%%%%%%%%%%%%%
    
    if(!silent){
      setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    
    ##%%%%%%%%%%%%%%
    ### UPDATE PARAMETERS
    ##%%%%%%%%%%%%%%
    ## extract BLUPs and BLUEs
    ulist <-list()
    for(k in 1:length(ind1g)){
      ulist[[k]] <- as.matrix(ge@x[ind1g[k]:ind2g[k]])
    }
    b <- as.matrix(ge@x[1:nx])
    Xb <- (X%*%b)
    
    Zulist <- list()
    Welist <- list()
    for(k in 1:length(ind1g)){
      Zulist[[k]] <- ZETA[[k]]$Z %*% ulist[[k]]
      Welist[[k]] <- Zulist[[k]]/varcomz[k]
    }
    zu <- do.call(cbind,Zulist); zu <- rowSums(zu)
    ##############################
    ## new estimate for error variance
    ## y'y - new.y'y / (N - nx)
    ##############################
    now <- 0 
    
    for(f in 1:length(ind1g)){ # y'Xb + y'Zu
      now <- now + crossprod(as.matrix(Zulist[[f]]),y)
    }
    b <- as.matrix(ge@x[ind1[1]:ind2[1]])
    now <- now + crossprod(X%*%b,y)
    
    varcomr[1] <- ( (t(y)%*%as.matrix(Rsei)%*%y) - now ) / (length(y)-nx)
    
    ##############################
    ## new estimates for variance components
    ##############################
    for(k in 1:length(ind1g)){ # adjust var comps except ERROR VARIANCE
      Kinv <- Ki[[k]]
      nk <- nrow(Kinv)
      rox <- ind1g[k]; cox <- ind2g[k]
      varcomz[k] <- ((t(ulist[[k]]) %*% Kinv  %*% ulist[[k]] ) + (as.numeric(varcomr[1])*sum(diag(Kinv%*%CMi[rox:cox,rox:cox]))))/nk 
    }
    
    #print(c(varcomz,varcomr))
    varcom <- c(varcomz,varcomr)
    zeros <- which(varcom <= 0)
    
    if(length(zeros)>0 & constraint){
      varcom[zeros] <- varcomr[1] * (1.011929e-07^2)
    }
    ## move back to the z and r separation
    names(varcom) <- c(rep("z",nz),rep("r",nr))
    varcomz <- varcom[which(names(varcom)=="z")]
    varcomr <- varcom[which(names(varcom)=="r")]
    #y <- Xb+zu
    varso <- cbind(varso,as.matrix(varcom))# just to keep track
  } ################# ======================== END OF ALGORITHM =======================
  ################### =================================================================
  ################### =================================================================
  ################### =================================================================
  e <- y - (Xb+zu)
  
  
  We <- (e)/varcomr[1]
  B <- cbind(do.call(cbind,Welist),We)
  left <- t(XZbind) %*% Rsei %*% B #left of the M matrix
  top <- t(left) # top of M
  corner <- t(B)%*%Rsei%*%B # left corner of M
  M <- rbind(cbind(corner,top),cbind(left,CM))
  vb <- try(chol(as(M, Class="sparseMatrix")), silent = TRUE)
  if(class(vb)=="try-error"){
    aii <- diag(nz+nr)
  }else{
    L <- t(vb); #L[1:4,1:4]
    ai <- L[1:(nz+nr),1:(nz+nr)]^2;ai
    ai <- as.matrix(ai)#*100#(var(y,na.rm = TRUE)) #originally should not be multiplied
    ai[upper.tri(ai)] <- t(ai)[upper.tri(ai)]
    aii <- solve(ai)
  }
  
  
  ## plot likelihood
  if(draw){
    layout(matrix(1:2,2,1))
    plot(y=ll.stored[-1], x=1:length(ll.stored[-1]), type="l", main="Log-likelihood", xlab="iteration",ylab="Log-likelihood")
    ## plot varcomp
    for(l in 1:dim(varso)[1]){
      #lb <- min(unlist(varso),na.rm = TRUE)
      ub <- max(unlist(varso),na.rm = TRUE)
      if(l==1){plot(varso[1,],ylim=c(0,ub),type="l", main="MME-EM results",ylab="varcomp", xlab="iteration")}else{
        lines(varso[l,],x=1:dim(varso)[2], col=l)
      }
    }
  }
  ## ======================== ##
  ## ======================== ##
  ## PROVIDE EXTRA PARAMETERS
  ## ======================== ##
  ## ======================== ##
  
  ## variance components
  out1 <- as.matrix(varcom); rownames(out1) <- varosss; colnames(out1) <- "component"
  ## inverse of the phenotypic variance (ZKZ'+R)-
  Vinv <- XZbind %*% CMi %*% t(XZbind)
  ## BLUPs
  names(ulist) <- varosssZ
  for(f in 1:length(ZETA)){
    rownames(ulist[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  ## VarBLUPs
  Var.u <- list()
  for(f in 1:length(ind1g)){
    Var.u[[f]] <- diag(CMi[ind1g[f]:ind2g[f],ind1g[f]:ind2g[f]])
    names(Var.u[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  ## betas
  rownames(b) <- colnames(X)
  ## var.betas
  xvxi <- CMi[ind1[1]:ind2[1],ind1[1]:ind2[1]]
  rownames(xvxi) <- colnames(xvxi) <- colnames(X)
  ## cond. residuals
  #e
  ## residuals
  ee <- y - (Xb)
  ## log likelihood
  ll <- as.vector(ll)
  ## AIC and BIC
  AIC = as.vector((-2 * ll ) + ( 2 * nx))
  BIC = as.vector((-2 * ll ) + ( log(length(y)) * nx))
  ## fish.inv
  # aii
  ## Ksp
  Ksp <- do.call(adiag1,lapply(ZETA,function(x){x$K}))
  ## fitted values only for non-missing data
  fit0 <- Xb + zu
  
  ## ======================== ##
  ## ======================== ##
  ## 2. PROVIDE EXTRA PARAMETERS
  ## using original data
  ## ======================== ##
  ## ======================== ##
  Xor.b <- X.or%*%b 
  Zor.u <-  lapply(as.list(1:nz),function(x,z,u){z[[x]]$Z %*% ulist[[x]]},z=ZETA.or,u=ulist)
  Zor.u2 <- do.call(cbind,Zor.u); Zor.u2 <- rowSums(Zor.u2)
  fit1 <- Xor.b + Zor.u2 # fitted values
  
  if(!silent){
    setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  }
  
  ### let user knows the constraints
  out1 <- as.data.frame(out1)
  out1[,"constraint"] <- "Positive"
  if(length(zeros)>0){
    out1[zeros,"constraint"] <- "Boundary"
    out1[zeros,1] <- 0
  }
  
  sigma.scaled <- out1[,1]/var(y,na.rm=TRUE)
  
  res <- list(var.comp=out1, V.inv = Vinv, u.hat=ulist, Var.u.hat=Var.u, 
              #PEV.u.hat=PEV.u, 
              beta.hat=b, Var.beta.hat=xvxi, residuals=ee, cond.residuals=e,
              LL=ll, sigma.scaled=sigma.scaled,
              AIC=AIC, BIC=BIC, fish.inv=aii,fitted.y.good=fit0, 
              X=X.or, Z=Zbind, K=Ksp, ZETA=ZETA,
              ## go back to original data
              fitted.y=fit1, fitted.u=Zor.u2, 
              forced=forced, convergence=convergence,
              CM=CM, CMi=CMi)
  
  layout(matrix(1,1,1))
  return(res)
}