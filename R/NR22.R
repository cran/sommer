NR22 <- function(y, X=NULL, ZETA=NULL, R=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=15, 
                constraint=TRUE, init=NULL, sherman=FALSE, che=TRUE, MTG2=FALSE, Fishers=FALSE, 
                gss=TRUE, forced=NULL, identity=TRUE, kernel=NULL, start=NULL, taper=NULL,
                verbose=0, gamVals=NULL, maxcyc=15, tol=1e-4){
  
  ## when two matrices are passed to regress this is also called
  ## to evaluate the REML at certain values of gamma and find a
  ## good place to start the regress algorithm
  SWsolve <- function(S,K,D,Dinv=NULL,b) {
    ## solve(a,b) where a has the form SKS' + D using the Sherman Morrison Woodbury identities where D is an easily inverted matrix
    
    if(is.matrix(K) & is.matrix(D) & !is.null(Dinv)) {
      ## Case 1 - all are matrices and D is already inverted
      tSDi <- crossprod(S,Dinv)
      Kinv <- solve(K)
      ret <- solve(Kinv + tSDi %*% S, tSDi)
      ret <- Dinv - crossprod(tSDi,ret)
      if(!missing(b)) ret <- ret %*% b
      return(ret)
    }
    
    if(is.numeric(K) & !is.null(Dinv)) {
      ## Case 2 - K is a number, ie first variance component is block random effect, random effects are iid
      tSDi <- crossprod(S,Dinv)
      ret <- solve(1/K * diag(ncol(S)) + tSDi %*% S, tSDi)
      ret <- Dinv - crossprod(tSDi,ret)
      if(!missing(b)) ret <- ret %*% b
      return(ret)
    }
    
    if(is.numeric(D) & is.matrix(K)) {
      ## diagonal D with a vector of entries along diagonal supplied
      ret <- 1/D * diag(nrow(S)) - 1/D^2 * S %*% solve(solve(K) + 1/D * crossprod(S),t(S))
      if(!missing(b)) ret <- ret %*% b
      return(ret)
    }
    
    if(is.numeric(K) & is.numeric(D)) {
      ret <- 1/D * diag(nrow(S)) - 1/D^2 * S %*% solve(1/K * diag(ncol(S)) + 1/D * crossprod(S),t(S))
      if(!missing(b)) ret <- ret %*% b
      return(ret)
    }
  }
  
  
  SWsolve2 <- function(Zlist,clist,b) {
    ## implementation of solve(Sigma,b) where Sigma is of the form sum [ clist[i] tcrossprod(Zlist[[ii]])) ] + diagonal with constant values along diagonal
    ## Invert a matrix of the form sum [ clist[i] tcrossprod(Zlist[[ii]])) ] using Sherman Morrison Woodbury Identities
    if(length(Zlist)!=(length(clist)-1)) stop()
    k <- length(Zlist)
    D <- clist[1] * tcrossprod(Zlist[[1]])
    diag(D) <- diag(D) + clist[k+1]
    Dinv <- SWsolve(Zlist[[1]],clist[1],clist[k+1])
    if(k==1) {
      if(!missing(b)) Dinv <- Dinv %*% b
      return(Dinv)
    }
    for(ii in 2:k) {
      Dinv <- SWsolve(Zlist[[ii]],clist[ii],D,Dinv)
      D <- D + clist[ii]*tcrossprod(Zlist[[ii]])
    }
    if(!missing(b)) Dinv <- Dinv %*% b
    return(Dinv)
  }
  
  reml <- function(lambda, y, X, V0, V1,verbose=0){
    
    if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)
      
      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
    V0 <- V0[isNA==F,isNA==F]
    V1 <- V1[isNA==F,isNA==F]
    X <- X[isNA==F,]
    X <- as.matrix(X)
    
    qr <- qr(X)
    ##print(qr$rank)
    n <- dim(as.matrix(y))[1]
    In <- diag(1,n)
    
    X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
    llik <- rep(0, length(lambda))
    if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]
    
    n <- dim(X)[1]
    if(missing(V0)) V0 <- diag(rep(1, n), n, n)
    rank <- n - qr$rank
    ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")
    for(i in 1:length(lambda))
    {
      if(verbose>=2) cat(lambda[i],"\n")
      
      Sigma <- (1-lambda[i])*V0 + lambda[i] * V1
      ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
      ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
      ##W <- chol2inv(cholesky)
      ##WX <- W %*% X
      WX <- solve(Sigma,cbind(X,In))
      W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      WX <- WX[,1:dim(X)[2]]
      XtWX <- t(X)%*%WX
      WQ <- W - WX%*%solve(XtWX,t(WX))
      rss <- t(y) %*% WQ %*% y
      logdetrss <- sum(log(eigen(rss)$values[1:q]))
      eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
      ldet <- sum(log(eVals))
      llik[i] <- Re(ldet*q/2 - rank*logdetrss/2)
    }
    imax <- sort.list(-llik)[1]
    lambdamax <- lambda[imax]
    curv <- 0
    if(imax > 1 && imax < length(lambda)){
      delta <- (lambda[imax+1] - lambda[imax-1])/2
      slope <-  (llik[imax+1] - llik[imax-1])/2
      curv <- llik[imax-1] -2*llik[imax] + llik[imax+1]
      lambdamax <- lambdamax - slope/curv * delta
      curv <- -curv/delta^2
    }
    lamMax <- lambdamax
    Sigma <- (1-lamMax)*V0 + lamMax * V1
    ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
    ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
    ##W <- chol2inv(cholesky)
    ##WX <- W %*% X
    WX <- solve(Sigma,cbind(X,In))
    W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
    WX <- WX[,1:dim(X)[2]]
    XtWX <- t(X)%*%WX
    FItWX <- solve(XtWX,t(WX))
    WQ <- W - WX%*%FItWX
    rss <- t(y) %*% WQ %*% y
    beta <- FItWX %*% y
    
    list(llik=as.numeric(llik),rms=rss/rank, beta=beta, gamma=lambda, gamMax=lambdamax,W=W)
  }
  
  remlOptimize <- function(y, X, V0, V1,verbose=0,...){
    
    if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)
      
      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
    V0 <- V0[isNA==F,isNA==F]
    V1 <- V1[isNA==F,isNA==F]
    X <- X[isNA==F,]
    X <- as.matrix(X)
    
    qr <- qr(X)
    ##print(qr$rank)
    n <- dim(as.matrix(y))[1]
    In <- diag(1,n)
    
    X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
    if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]
    
    n <- dim(X)[1]
    if(missing(V0)) V0 <- diag(rep(1, n), n, n)
    rank <- n - qr$rank
    ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")
    
    f <- function(lambda,verbose=verbose) {
      if(verbose>=2) cat(lambda,"\n")
      Sigma <- (1-lambda)*V0 + lambda * V1
      ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
      ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
      ##W <- chol2inv(cholesky)
      ##WX <- W %*% X
      WX <- solve(Sigma,cbind(X,In))
      W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      WX <- WX[,1:dim(X)[2]]
      XtWX <- t(X)%*%WX
      WQ <- W - WX%*%solve(XtWX,t(WX))
      rss <- t(y) %*% WQ %*% y
      logdetrss <- sum(log(eigen(rss)$values[1:q]))
      eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
      ldet <- sum(log(eVals))
      llik <- Re(ldet*q/2 - rank*logdetrss/2)
      llik
    }
    
    res <- optimize(f,interval=c(0,1),maximum=TRUE,verbose=verbose,...)
    lamMax <- res$maximum
    llikMax <- res$objective
    
    Sigma <- (1-lamMax)*V0 + lamMax * V1
    ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
    ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
    ##W <- chol2inv(cholesky)
    ##WX <- W %*% X
    WX <- solve(Sigma,cbind(X,In))
    W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
    WX <- WX[,1:dim(X)[2]]
    XtWX <- t(X)%*%WX
    FItWX <- solve(XtWX,t(WX))
    WQ <- W - WX%*%FItWX
    rss <- t(y) %*% WQ %*% y
    beta <- FItWX %*% y
    
    list(llik=as.numeric(llikMax),rms=rss/rank, beta=beta, gamMax=lamMax,W=W)
  }
  
  
  ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #ZETA.or <-ZETA
  if(is.null(names(ZETA))){
    varosss <- c(paste("u.",1:length(ZETA), sep=""))
  }else{
    varosss <- as.character(names(ZETA))
  }
  
  # model has fixed and random effects
  if(!is.null(X) & !is.null(ZETA)){ 
    KKK <- list()
    for(h in 1:(length(ZETA))){
      if(is.square.matrix(ZETA[[h]][[1]])){
        if(is.diagonal.matrix(ZETA[[h]][[1]])){
          KKK[[h]] <- (ZETA[[h]][[2]])
        }else{
          KKK[[h]] <- tcrossprod(ZETA[[h]][[1]], ZETA[[h]][[1]] %*% (ZETA[[h]][[2]]) ) 
        }
      }else{
        KKK[[h]] <- tcrossprod(ZETA[[h]][[1]], ZETA[[h]][[1]] %*% (ZETA[[h]][[2]]) ) 
      }
    }
    formula = as.formula(y~X) 
    Vformula = as.formula(paste("~ ",paste(paste("KKK[[",1:length(ZETA),"]]",sep=""), collapse="+")))
  }
  else if(is.null(X) & !is.null(ZETA)){
    KKK <- list()
    for(h in 1:(length(ZETA))){
      if(is.square.matrix(ZETA[[h]][[1]])){
        if(is.diagonal.matrix(ZETA[[h]][[1]])){
          KKK[[h]] <- (ZETA[[h]][[2]])
        }else{
          KKK[[h]] <- tcrossprod(ZETA[[h]][[1]], ZETA[[h]][[1]] %*% (ZETA[[h]][[2]]) ) 
        }
      }else{
        KKK[[h]] <- tcrossprod(ZETA[[h]][[1]], ZETA[[h]][[1]] %*% (ZETA[[h]][[2]]) ) 
      }
    }
    formula = as.formula(y~1) 
    Vformula = as.formula(paste("~ ",paste(paste("KKK[[",1:length(ZETA),"]]",sep=""), collapse="+")))
  }
  
  if(!is.null(X)){
    x.or <- X
  }else{
    x.or <- matrix(rep(1,length(y)))
  }
  
  
  ## Vformula can just be something like ~ V0 + V1
  ## or leave it out or Vformula=NULL
  ## assume its in the form ~ V1 + V2 + ... + Vn or missing or Vformula=NULL
  ## for random effects and random interactions for factors A and B include
  ## ~ A + B + I(A:B)
  
  if(verbose>9) cat("Extracting objects from call\n")
  #if(missing(data)) 
  data <- environment(formula)
  mf <- model.frame(formula,data=data,na.action=na.pass)
  mf <- eval(mf,parent.frame())
  y <- model.response(mf)
  
  model <- list()
  model <- c(model,mf)
  
  if(missing(Vformula)) Vformula <- NULL
  
  ## Find missing values in fixed part of the model :: Change Aui Mar 1 2012
  isNA <-  apply(is.na(mf), 1, any)
  
  if(!is.null(Vformula))
  {
    V <- model.frame(Vformula,data=data,na.action=na.pass)
    V <- eval(V, parent.frame())
    # find missings in random part of the model Aui Mar 1 2012
    mfr <- is.na(V)
    if(ncol(mfr) == 1){
      isNA <- isNA | mfr
    } else {
      isNA <- isNA | apply(mfr[,!apply(mfr,2,all)], 1, any) # use only columns of the matrix with some nonmissing values
    }
    rm(mfr)
    Vcoef.names <- names(V)
    V <- as.list(V)
    k <- length(V)
  } else {
    V <- NULL
    k <- 0
    Vcoef.names=NULL
  }
  
  if(ncol(mf)==1) mf <- cbind(mf,1)
  X <- model.matrix(formula, mf[!isNA,]) # Aui Mar 1 2012 account for missings in the random part
  
  y <- y[!isNA]
  n <- length(y)
  Xcolnames <- dimnames(X)[[2]]
  if(is.null(Xcolnames)) {
    Xcolnames <- paste("X.column",c(1:dim(as.matrix(X))[2]),sep="")
  }
  
  X <- matrix(X, n, length(X)/n)
  qr <- qr(X)
  rankQ <- n-qr$rank
  if(qr$rank) {
    X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
    Xcolnames <- Xcolnames[qr$pivot[1:qr$rank]]
  } else {
    cat("\nERROR: X has rank 0\n\n")
  }
  
  if(verbose>9) cat("Setting up kernel\n")
  
  if(missing(kernel)){
    K <- X
    colnames(K) <- Xcolnames
    reml <- TRUE
    kernel <- NULL
  } else {
    if(length(kernel)==1 && kernel>0){
      K <- matrix(rep(1, n), n, 1)
      colnames(K) <- c("1")
    }
    if(length(kernel)==1 && kernel<=0) {
      K <- Kcolnames <- NULL
      KX <- X
      rankQK <- n
    }
    if(length(kernel) > 1) {
      ##K is a matrix I hope :: Change Aui Mar 1 2012
      if(is.matrix(kernel)) {
        K <- kernel[!isNA,]
      } else {
        K <- model.frame(kernel, data=data, na.action=na.pass)
        K <- eval(K, parent.frame())
        if(ncol(K) == 1){
          dimNamesK <- dimnames(K)
          K <- K[!isNA, ]
          dimNamesK[[1]] <- dimNamesK[[1]][!isNA]
          K <- data.frame(V1 = K)
          dimnames(K) <- dimNamesK
        } else {
          K <- K[!isNA, ]
        }
        K <- model.matrix(kernel, K)
      }
    }
    reml <- FALSE
  }
  
  if(!is.null(K)){
    Kcolnames <- colnames(K)
    qr <- qr(K)
    rankQK <- n - qr$rank
    if(qr$rank == 0) K <- NULL else {
      K <- matrix(K[, qr$pivot[1:qr$rank]],n,qr$rank)
      Kcolnames <- Kcolnames[qr$pivot[1:qr$rank]]
      KX <- cbind(K, X) # Spanning K + X: Oct 12 2011
      qr <- qr(KX)
      KX <- matrix(KX[, qr$pivot[1:qr$rank]],n,qr$rank) # basis of K+X
    }
  }
  
  if(missing(maxcyc)) maxcyc <- 50
  if(missing(tol)) tol <- 1e-4
  delta <- 1
  
  if(verbose>9) cat("Removing parts of random effects corresponding to missing values\n")
  ## remove missing values
  for(i in 1:k){
    if(is.matrix(V[[i]]))
    {
      V[[i]] <- V[[i]][!isNA, !isNA]
    }
    if(is.factor(V[[i]]))
    {
      V[[i]] <- V[[i]][!isNA]
    }
  }
  
  In <- diag(rep(1,n),n,n)
  
  if(identity) {
    V[[k+1]] <- as.factor(1:n)
    names(V)[k+1] <- "In"
    k <- k+1
    
    Vcoef.names <- c(Vcoef.names,"In")
    Vformula <- as.character(Vformula)
    Vformula[1] <- "~"
    Vformula[2] <- paste(Vformula[2],"+In")
    Vformula <- as.formula(Vformula)
  }
  
  model <- c(model,V)
  model$formula <- formula
  model$Vformula <- Vformula
  
  ## specify which parameters are positive and which are negative
  ## pos = c(1,1,0) means first two parameters are positive, third is either
  
  #if(!missing(pos)) pos <- as.logical(pos)
  pos <- rep(FALSE,k)
  #if(length(pos) < k) cat("Warning: argument pos is only partially specified; additional terms (n=",k-length(pos),") set to FALSE internally.\n",sep="")
  #pos <- c(pos,rep(FALSE,k))
  #pos <- pos[1:k]
  ## Sherman Morrison Woodbury identities for matrix inverses can be brought to bear here
  if(verbose>9) cat("Checking if we can apply the Sherman Morrison Woodbury identites for matrix inversion\n")
  if (all(sapply(V, is.factor)) & k>2 ) {  # Contribution by Hans Jurgen Auinger
    SWsolveINDICATOR <- TRUE
  } else SWsolveINDICATOR <- FALSE
  #print(SWsolveINDICATOR)
  Z <- list()
  for (i in 1:length(V)) {
    if (is.factor(V[[i]])) {
      Vi <- model.matrix(~V[[i]] - 1)
      colnames(Vi) <- levels(V[[i]])
      Z[[i]] <- Vi
      V[[i]] <- tcrossprod(Vi)
    } else{
      Z[[i]] <- V[[i]]
    }
  }
  names(Z) <- names(V)
  
  ## So V is always a list of variance coavriance matrices, Z contains
  ## the model matrices of factors when we need to invoke the Sherman
  ## Woodbury identities
  
  ## Expected Fisher Information
  A <- matrix(rep(0, k^2), k, k)
  entries <- expand.grid(1:k,1:k)
  
  x <- rep(0,k)
  sigma <- c(1,rep(0, k-1))
  
  stats <- rep(0, 0)
  
  ## START ALGORITHM
  #print(taper)
  # taper is the steps to take 0-1
  #print(taper)
  if(missing(taper)){
    taper <- rep(0.9, maxcyc)
    if(missing(start) && k>1) taper[1:2] <- c(0.5, 0.7)
  } else {
    taper <- pmin(abs(taper), 1)
    if((l <- length(taper)) < maxcyc) taper <- c(taper, rep(taper[l], maxcyc-l))
  }
  
  if(!is.null(start)) {
    ## pad start with zeros if required
    start <- c(start, rep(1,k))
    start <- start[1:k]
  }
  # starting values for the random effects
  if(k>2 && is.null(start)) start <- rep(var(y,na.rm=TRUE),k)
  if(k==1 && is.null(start)) start <- var(y,na.rm=TRUE)
  
  # skip this
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(is.null(start) && k==2) {
    if(missing(gamVals)) {
      gamVals <- seq(0.01,0.02,length=3)^2
      gamVals <- sort(c(gamVals,seq(0.1,0.9,length=3),1-gamVals))
      gamVals <- 0.5
    }
    if(length(gamVals)>1) {
      if(verbose>=1) cat("Evaluating the llik at gamma = \n")
      if(verbose>=1) cat(gamVals)
      if(verbose>=1) cat("\n")
      reg.obj <- reml(gamVals,y,X,V[[1]],V[[2]],verbose=verbose)
      llik <- reg.obj$llik
      llik <- as.double(llik)
      if(verbose>=2) cat(llik,"\n")
      gam <- gamVals[llik==max(llik)]
      gam <- gam[1]
      if(verbose>=2) cat("MLE is near",gam,"and llik =",max(llik),"there\n")
    }
    if(length(gamVals)==1) {
      ## go straight to the Newton Raphson at gamVals
      gam <- gamVals[1]
      reg.obj <- list(rms=var(y))
    }
    start <- c(1-gam,gam)*reg.obj$rms
    ## it tends to take huge steps when starting at gam=0.9999
    if(gam==0.9999) {
      taper[1] <- taper[1]/100
      maxcyc <- maxcyc*10
    }
    if(verbose>=1) cat(c("start algorithm at",round(start,4),"\n"))
  }
  
  if(is.null(start) & k>2) {
    ## Never gets here by default - but this could be implemented,
    ## though it does add on a few extra iterations at the
    ## start.... not necessary in basic examples
    
    LLvals <- NULL
    ## equal weights
    V2 <- V[[2]]
    for(ii in 3:k) V2 <- V2 + V[[ii]]
    LLvals <- c(LLvals,reml(0.5,y,X,V[[1]],V2)$llik)
    ## Most at one end
    V2 <- V[[1]] + V2 ## total of all Vs
    for(ii in 1:k) {
      V2 <- V2 - V[[ii]]
      LLvals <- c(LLvals,reml(0.75,y,X,V2,V[[ii]])$llik)
    }
    best <- which.max(LLvals)
    if(verbose) {
      cat("Checking starting points\n")
      cat("llik values of", LLvals, "\n")
    }
    if(best==1) {
      start <- rep(var(y,na.rm=TRUE),k)
    } else {
      start <- rep(0.25,k)
      start[best] <- 0.75
    }
  }
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(is.null(init)){
    sigma <- coef <- start
  }else{
    sigma <- coef <- init
  }
  ## reparameterise so everything will get into the correct spot after exp
  coef[pos] <- log(sigma[pos])
  coef[!pos] <- sigma[!pos]
  
  ## Set the memory requirements beforehand
  T <- vector("list", length=k)
  for(ii in 1:k) T[[ii]] <- matrix(NA,n,n)
  
  logL2.stored <- numeric()
  record <- as.matrix(sigma)#matrix(NA,nrow=length(sigma), ncol=)
  
  ##################
  ## progress bar
  if(!silent){
    count <- 0
    tot <- 15
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  ################## ==========================================================================
  ################## ==========================================================================
  ################## ==========================================================================
  ################## ======================= START VAR.COMP ESTIMATION ========================
  for(cycle in 1:maxcyc){
    
    if(!silent){
      count <- count + 1
    }
    ## Limit how far out we go on the logarithmic scale
    ind <- which(pos)
    #print(ind)
    if(length(ind)) {
      coef[ind] <- pmin(coef[ind],20)
      coef[ind] <- pmax(coef[ind],-20) ## so on regular scale everything is between exp(-20) and exp(20)
      sigma[ind] <- exp(coef[ind])
    }
    
    if(verbose) {
      cat(cycle, "sigma =",sigma)
      ##cat(sigma)
    }
    #print(SWsolveINDICATOR)
    if(!SWsolveINDICATOR) {
      #print("yeah")
      Sigma <- 0
      #print(V[[1]][1:5,1:5])
      ## can we get rid of this loop? multiplies each var-cov for its variance component
      ## additionally sums up each var-cov V = ZKZ + ZKZ + ... + ZRZ
      for(i in 1:k) {
        #print(V[[i]][1:3,1:5])
        #print(sigma[i])
        Sigma <- Sigma + V[[i]]*sigma[i]
      }
      #print(Sigma[1:4,1:4])
      #print(k)
      #print(sigma)
      #for(i in 1:k) print(V[[i]][1:5,1:5])
      # solves V to get V.inverse using In (identity matrix) as response
      W <- solve(Sigma,In)
      #print(W[1:3,1:3])
    } else { # if no Vformula
      W <- SWsolve2(Z[1:(k-1)],sigma)
    }
    # K is the X matrix for fixed effects
    if(is.null(K)){
      WQK <- W
    }else{ # our case, use sherman
      WK <- W %*% K # V- X
      WQK <- W - WK %*% solve(t(K)%*%WK, t(WK)) # P = V-  -  V-X [X' V- X]-1 XV-
    }
    ## WQK and WQX
    if(reml){ # REML(default)
      WQX <- WQK # WQK and WKX are the same
    }else{ # ML(not default)
      WX <- W %*% KX        # including the kernel (Oct 12 2011)
      WQX <- W - WX %*% solve(t(KX)%*%WX, t(WX))
    }
    
    # y' P y
    rss <- as.numeric(t(y) %*% WQX %*% y)
    
    ##if(verbose>9) cat("Sigma[1:5]",Sigma[1:5],"\n")
    ##if(verbose>9) cat("RSS",rss,"WQX[1:5]",WQX[1:5],"\n")
    sigma <- sigma * rss/rankQK
    coef[!pos] <- sigma[!pos]
    coef[pos] <- log(sigma[pos])
    WQK <- WQK * rankQK/rss # P * [r(X) / yPy]
    WQX <- WQX * rankQK/rss
    rss <- rankQK ## looks bad but the rss is absorbed into WQK so the rss term comes out of eig below
    
    # this is the way to check if the 
    eig <- sort(eigen(WQK,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankQK]
    
    if(any(eig < 0)){
      WQK <- WQK + (tol - min(eig))*diag(dim(WQK)[1])
      eig <- eig + tol - min(eig)
    }
    ldet <- sum(log(eig))
    llik <- ldet/2 - rss/2
    if(cycle == 1) llik0 <- llik
    delta.llik <- llik - llik0
    llik0 <- llik
    
    ## From Jean-Luc Jannick to fix excess carriage returns
    if (verbose){
      if (reml) cat(" resid llik =", llik, "\n") else cat(" llik =", llik, "\n")
      cat(cycle, "adjusted sigma =", sigma)
      if (cycle > 1){
        if (reml) cat(" delta.llik =", delta.llik, "\n") else cat(" delta.llik =", delta.llik, "\n")
      } else cat("\n")
    }
    
    ## now the fun starts, derivative and expected fisher info
    ## the 0.5 multiple is ignored, it is in both and they cancel
    
    ##T <- list(NULL)
    x <- NULL
    
    ## derivatives are now D[[i]] = var.components[i]*V[[i]]
    var.components <- rep(1,k)
    ind <- which(pos)
    if(length(ind)) var.components[ind] <- sigma[ind]
    
    ## Slow part - order k n-squared P Vi
    #print(identity)
    if(!SWsolveINDICATOR) {
      if(identity) {
        T[[k]] <- WQK
        if(k>1) {
          for(ii in (k-1):1) T[[ii]] <- WQK %*% V[[ii]] ## # P ZKZ' #... V.inv - V.inv X [X' V.inv X]-1 * ZKZi
        }
      } else {
        for(ii in 1:k) T[[ii]] <- WQK %*% V[[ii]]
      }
    } else {
      if(identity) {
        T[[k]] <- WQK
        if(k>1) {
          for(ii in (k-1):1) T[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]]) # PZZ'
        }
      } else {
        for(ii in 1:k) T[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
      }
    }
    ## obtain first derivatives
    # Vi = dV/ds
    # V.inv - V.inv X [X' V.inv X]-1, adjust +- variance components
    ## y' P Vi P y - tr(P Vi)
    x <- sapply(T,function(x) as.numeric(t(y) %*% x %*% WQX %*% y - sum(diag(x))))
    ## theta(k) * 1st.deriv are scalar values
    x <- x * var.components
    
    
    ## See nested for loops commented out below - 
    ## evaluating the Expected Fisher Information, A
    ## [theta(i) * 1st.deriv(i)] * [theta(j) * 1st.deriv(j)]  * sigma(i) * sigma(j)
    #print(var.components)
    ff <- function(x) sum(T[[x[1]]] * t(T[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
    aa <- apply(entries,1,ff) # matrix of combinations of var.comp 
    ## 1 1
    ## 2 1
    ## 1 2
    ## 2 2
    # A is the Fishers information matrix
    A[as.matrix(entries)] <- aa
    
    
    stats <- c(stats, llik, sigma[1:k], x[1:k])
    if(verbose>=9) {
      ##cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
    }
    
    logL2.stored[cycle] <- llik
    record <- cbind(record,as.matrix(sigma))
    A.svd <- ginv(A)
    x <- A.svd %*% x
    
    if(qr(A)$rank < k){
      if(cycle==1) {
        if(verbose) {
          cat("Warning: Non identifiable dispersion model\n")
          ##print(round(A,6))
          cat(sigma)
          cat("\n")
        }
      }
    }
    
    ## end of newton-raphson step
    ## x is  -l'(sigma)/l''(sigma)
    ## hence we add instead of subtract
    ## taper controls the proportion of each step to take
    ## for some reason 0.5 works very well
    
    ##if(all(pos==1)) x <- sign(x) * pmin(abs(x),5) ## limit maximum shift we can take in one step to 5 units on log scale
    coef <- coef + taper[cycle] * x
    sigma[!pos] <- coef[!pos]
    sigma[pos] <- exp(coef[pos])
    
    if(draw){# draw
      ylim <- max(unlist(record), na.rm=TRUE)
      my.palette <- brewer.pal(7,"Accent")
      layout(matrix(1:2,2,1))
      plot(logL2.stored,type="l", main="logLikelihood", col=my.palette[7],lwd=3, las=2, xaxt="n", ylab="logLikelihood value", xlab="Iterations processed", cex.axis=0.5) 
      axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
      legend("bottomright", legend = round(llik,3), bty="n", cex=0.7)
      plot(record[1,],ylim=c(0,ylim),type="l", las=2, xaxt="n",main="Newton-Raphson algorithm results", col=my.palette[1],lwd=3, ylab="Value of the variance component", xlab="Iterations processed", cex.axis=0.5) 
      axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
      for(t in 1:(dim(record)[1])){
        lines(record[t,],col=my.palette[t],lwd=3)
      } 
      ww <- dim(record)[1]
      lege <- list()
      #lege2 <- list()
      for(ka in 1:length(sigma)){
        if(ka == length(sigma)){
          lege[[ka]] <- paste("Var(e):",round(record[ka,cycle+1],4), sep="")
          #lege2[[k]] <- paste("Var(e):")
        }else{
          lege[[ka]] <- paste("Var(",varosss[ka],"):",round(record[ka,cycle+1],4), sep="")
          #lege2[[k]] <- paste("Var(u",k,"):",sep="")
        }
      }
      legend("topright",bty="n", cex=0.7, col=my.palette, lty=1, lwd=3, legend=unlist(lege))
      
    }########end draw
    
    if(!silent){
      setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    
    ###### controls
    #vv <- which(sigma < 0)
    #if(length(vv)>0){
    #  sigma[vv] <- abs(sigma[vv])
    #}
    
    
    ## check the change in llik is small
    if(cycle > 1 & abs(delta.llik) < tol*10) {
      if(!silent){
        setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
      }
      break
    }
    if(max(abs(x)) < tol) {
      if(!silent){
        setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
      }
      break
    }
  }
  
  if(any(eig < 0)){
    if(!silent){
      cat("\nError: Sigma is not positive definite on contrasts: range(eig)=", range(eig), "\n")
    }
  }
  ################
  ############### ========================= FINISH VAR.COMP ESTIMATED =============== #############
  ############### ========================= ========================= ========================= 
  ############### ========================= ========================= ========================= 
  ############### ========================= ========================= ========================= 
  
  
  return(list(vars=sigma,W=W))
}
