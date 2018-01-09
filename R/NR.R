MNR <- function(Y,X=NULL,ZETA=NULL,R=NULL,init=NULL,iters=20,tolpar=1e-3,
                tolparinv=1e-6,draw=FALSE,silent=FALSE, constraint=TRUE, 
                EIGEND=FALSE, forced=NULL, IMP=FALSE, complete=TRUE, 
                check.model=FALSE, restrained=NULL, REML=TRUE, init.equal=TRUE){
  
  up.to.vec <- function(x){
    if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
    babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
    x2 <- x[babas]; names(x2) <- paste(rownames(x)[babas[,1]],rownames(x)[babas[,2]],sep=".")
    }else{x2 <- as.matrix(x)}
    return(x2)
  }
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  copying <- function(m) { # copy upper triangular in lower triangular
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  copying2 <- function(m) { # copy lower triangular in upper triangular
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  ############################
  if(check.model){ # if coming from mmer don't check
    if(is.list(ZETA)){
      if(is.list(ZETA[[1]])){ # if was provided as a two level list
        ZETA=ZETA
      }else{ # if was provided as a one level list
        ZETA=list(ZETA)
      }
    }else{
      #stop;
      stop(cat("\nThe random effects need to be provided in a list format, please see examples"), call. = FALSE)
    }
    ZETA <- lapply(ZETA, function(x){
      if(length(x) == 1){
        provided <- names(x)
        if(provided == "Z"){
          y <- list(Z=x[[1]],K=diag(dim(x[[1]])[2]))
        }else if(provided == "K"){
          y <- list(Z=diag(length(y)),K=x[[1]])
        }else{
          stop(call.=FALSE)
          cat("Names of matrices provided can only be 'Z' or 'K', the names you provided don't match the arguments required")
        }
      }else{y <- x}; 
      return(y)
    })
  }
  ###########################
  if(EIGEND){
    IMP <- TRUE
    cat("EIGEND feature activated. Eigen decomposition of K will be performed\n")
    Vmat.type <- "sparseMatrix" # if K will be sparse then V will be sparse
  }else{Vmat.type <- "dgeMatrix"} # if K will be dense V will be dense to be used in the algorithm
  Y <- as.matrix(Y); 
  if(is.null(colnames(Y))){colnames(Y) <- paste("T",1:ncol(Y),sep="")}
  traitnames <- colnames(Y)
  if(IMP){
    Y <- as.matrix(apply(Y,2,imputev))#
  }
  if(complete){# even IMP=FALSE imputes incomplete cases
    nona <- sort(unique(as.vector(unlist(apply(as.matrix(Y),2,function(x){which(!is.na(x))})))))
    Y[nona,] <- apply(as.matrix(Y[nona,]),2,imputev)
  }
  ###########################
  convergence <- FALSE
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### fix that some initially were not specified
  if(is.null(X)){
    X <- matrix(1,nrow=nrow(as.matrix(Y)))
    colnames(X) <- "(Intercept)"
  }
  if(is.null(R)){
    R <- list(units=diag(nrow(as.matrix(Y))))
  }
  fixnames <- colnames(X)
  Y.or <- Y
  X.or <- X # completly original same in dimensions
  #qr.Xor <- qr(X.or)
  
  ZETA.or <- ZETA
  R.or <- R
  if(is.null(names(ZETA))){
    varosss <- c(paste("u",1:length(ZETA), sep=""))
  }else{
    varosss <- c(names(ZETA))
  }; varosssZ <- varosss
  if(is.null(names(R))){
    varosss <- c(varosss,paste("Res",1:length(R),sep=""))
  }else{
    varosss <- c(varosss,names(R))
  }
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  ## if EIGEND (data is complete no missing data)
  if(EIGEND==TRUE){
    dias <- unlist(lapply(ZETA, function(x){is.diagonal.matrix(x$K)}))
    
    EIGENS <- eigen(ZETA[[1]]$K) # ONLY FOR THE 1ST K (only one covariance mat is allowed)
    Us <- EIGENS$vectors # extract eigen vectors U
    ###Usp <- as(do.call("adiag1", Us),Class="sparseMatrix") # U'G as diagonal
    Ds <- diag(EIGENS$values) # extract eigen values D
    ###Dsp <- as(do.call("adiag1", Ds),Class="sparseMatrix") # U'G as diagonal
    # transform ZETA
    ZETA[[1]]$K <- Ds
    tUs <- t(Us)
    tUsi <- solve(tUs) # inverse
    # controls
    if(length(which(!dias)) > 1){
      stop("Eigen decomposition only works for one dense relationship matrix.", call. = FALSE)
    }
    if(which(dias != TRUE) != 1){ # FALSE should be the 1st position
      stop("The random effect corresponding to the relationship matrix needs to be in the first position.", call. = FALSE)
    }
    if(dim(X)[1]!=dim(Us)[1]){
      stop("Eigen decomposition only works for square models (non-replicated experiments).", call. = FALSE)
    }
    # transformX
    X <- tUs %*% X
    # transform other Z's
    if(length(ZETA)>1){
      ZETA[-c(1)] <- lapply(ZETA[-c(1)], function(x){x$Z <- tUs%*%x$Z; return(x)})
    }
    # transform y
    Y <- tUs %*% as.matrix(Y)
  }
  Y.for.eigend <- Y
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### now reduce the original inputs
  good0 <- apply(as.matrix(Y),2,function(x){which(!is.na(x))})
  if(is.matrix(good0)){good<- sort(Reduce(intersect,split(good0, rep(1:ncol(good0), each = nrow(good0)))))
  }else{good<- sort(Reduce(intersect,good0))}
  if(!complete){nona<-good} # to store which where the indexes of the reduced dataset
  Y <- as.matrix(Y[good,]); colnames(Y) <- traitnames
  Y.red.noscale.ordim <- Y # Y reduced no scaled original dimensions
  Y <- scale(Y); colnames(Y) <- traitnames # Y reduced but scaled
  Y.red.scale.ordim <- Y
  X <- as.matrix(X[good,])
  X.noeigen <- as.matrix(X.or[good,])
  X.red.ordim <- X
  ZETA <- lapply(ZETA, function(x,good){x$Z <- x$Z[good,]; x$K<- x$K; return(x)}, good=good)
  R <- lapply(R, function(x,good){x <- x[good,good]; return(x)}, good=good)
  #### get sizes of reduced
  nz <- length(ZETA)
  nr <- length(R)
  #nx <- dim(X)[2]
  n <- nrow(Y)
  #k <- nz+nr
  #In <- as(diag(n),Class="sparseMatrix")
  qr <- qr(X)
  #ranx <- qr$rank # length(good)
  #rankQ <- n-qr$rank
  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  ## move everything to sparse matrices
  #print(str(ZETA))
  ZETA <- lapply(ZETA, function(x){ 
    x$Z <- as(x$Z,Class="sparseMatrix");
    x$K <- as(x$K,Class="sparseMatrix");
    return(x)
  })
  R <- lapply(R, function(x){as(x,Class="sparseMatrix")})
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  ## form c(ZKZ,R)
  ZKZ <- lapply(ZETA,function(x){
    if(is.square.matrix(as.matrix(x$Z))){
      if(is.diagonal.matrix(x$Z)){
        res <- (x$K)
      }else{
        if(length(x)==3){ # user gave Zt for example for unstructured models
          res <- tcrossprod(x$Zt, x$Z %*% (x$K) ) 
        }else{ # user only provided Z and K
          res <- tcrossprod(x$Z, x$Z %*% (x$K) )  
        }
      }
    }else{ # if z$Z is not square
      if(length(x)==3){ # user gave Zt for example for unstructured models
        res <- tcrossprod(x$Zt, x$Z %*% (x$K) ) 
      }else{ # user only provided Z and K
        res <- tcrossprod(x$Z, x$Z %*% (x$K) )  
      }
    }
    return(res)
  })
  ZKZ <- c(ZKZ,R)
  
  #no. of individuals originally
  if(IMP==FALSE){
    n <- dim(ZETA.or[[1]]$Z)[1]#no. of individuals in K matrix
  }else{n <- dim(ZETA[[1]]$Z)[1]} # not sure the else is needed ( i have modified too much)
  
  dimos <- dim(as.matrix(Y))
  ts <- dimos[2]
  inds <- dimos[1] #no. of individuals
  base.var <- var(as.matrix(Y.or),na.rm=TRUE) # original var
  if(EIGEND){base.var <- var(as.matrix(Y.for.eigend),na.rm=TRUE)} # for EIGEND we need the variance after transformation
  sc.var <- var(as.matrix(Y),na.rm=TRUE)# scaled var
  
  ###############
  # LINEARIZE
  ###############
  # impute phenotypes
  # decompose in a vector
  Y <- as.matrix(as.vector(Y)); #colnames(Y) <- traitnames #dim(Y) # Y scaled and linearized
  Y.red.noscale.lin <- as.matrix(as.vector(Y.red.noscale.ordim)) # Y reduced no scaled linearized
  tY <- t(Y); dim(tY) #get transpose
  X <- do.call("adiag1", rep(list(X), ts)); #dim(X) # X multivariate
  X.red.lindim <- X
  X.red.lin.noeigen <- do.call("adiag1", rep(list(X.noeigen), ts));# dim(X) # X multivariate
  tX <- t(X); #dim(tX)
  qr <- qr(X)
  rankX <- dim(X)[1]-qr$rank; #rankX
  
  ###############
  # INITIAL VAR VALUES
  ###############
  ## define variance-covariance components, each list is a var.comp
  if(is.null(init)){
    sigma <- rep(list(sc.var),nz+nr) #original variance values
    if(init.equal){ # if user prefers to use equal starting values (preferred-asreml)
      for(w in 1:length(sigma)){
        if(w <= nz){sigma[[w]] <- (sigma[[w]]*0 + 1)*0.10 + diag(0.05,ts)} # random
        if(w > nz){sigma[[w]] <- (sigma[[w]]*0 + 1)*0.04977728 + diag(0.02488864,ts)} #rcov
      }
    }
  }else{
    sigma <- init
  }
  if(!is.null(forced)){
    iters <- 1
    sigma <- forced
  }
  #print(sigma)
  names(sigma) <- varosss
  # decompose in a vector
  varos <- lapply(sigma,up.to.vec)
  sigma2 <- as.matrix(unlist(varos)) # current values of vc
  coef2 <- sigma2
  sigma3 <- sigma2
  pos <- rep(FALSE,length(sigma2))
  
  taper <- rep(0.9, iters) # weighting parameter
  taper[1:2] <- c(0.5, 0.7)
  k <- length(sigma2)
  llstore <- numeric()
  ###################
  # MAT FOR DERIVS
  ###################
  ## trait combinations
  traitm <- expand.grid(1:ts,1:ts)
  if(ts > 1){
    traitm <- (traitm[!duplicated(t(apply(traitm, 1, sort))),])[,c(2,1)]; colnames(traitm) <- c("t1","t2")
  }else{
    traitm <- as.matrix(cbind(traitm,traitm)[,1:2]) # when single trait we have issues
    colnames(traitm) <- c("t1","t2")
  }
  pos.mats <- list()
  ## for each trait-combo fill a dummy matrix for derivatives
  namos <- rownames(sigma[[1]]) #names of traits
  for(i in 1:dim(traitm)[1]){ 
    temp.mat <- matrix(0,ts,ts); 
    i1 <- traitm[i,1]
    i2 <- traitm[i,2]
    colnames(temp.mat) <- rownames(temp.mat) <- namos
    temp.mat[i1,i2] <- 1
    temp.mat[i2,i1] <- 1
    pos.mats[[i]] <- temp.mat
    names(pos.mats)[i] <- paste("Deriva",namos[i1],namos[i2], sep=".")
  }
  posmats.list.vc <- rep(list(pos.mats),nz+nr) # THIS DERIVATIVES ARE NEEDED FOR EVERY RANDOM EFFECT
  names(posmats.list.vc) <- varosss
  
  ####################
  ## mapper to know which terms un the linearized var.comp are diagonal
  diagss <- which(traitm[,1] == traitm[,2]) # this terms are diagonal or ti-ti
  offdiagss <- which(traitm[,1] != traitm[,2]) # this terms are diagonal or ti-ti
  diagss <- rep(list(diagss), nz+nr)
  offdiagss <- rep(list(offdiagss), nz+nr)
  all <- rep(list(traitm), nz+nr)
  for(u in 2:length(diagss)){
    diagss[[u]] <- diagss[[u]] + diagss[[u-1]][length(diagss[[u-1]])]
    offdiagss[[u]] <- offdiagss[[u]] + offdiagss[[u-1]][length(offdiagss[[u-1]])]
    all[[u]] <- all[[u]] + max(all[[u-1]])
  }
  names(diagss) <- varosss
  names(offdiagss) <- varosss
  diagss2 <- unlist(diagss)
  offdiagss2 <- unlist(offdiagss)
  all2 <- do.call(rbind,all); rownames(all2) <- NULL
  ####################
  var.comp.ret<- list()
  convergence <- FALSE
  if(!silent){
    count <- 0
    tot <- iters
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  ## ============================================================================
  ## ============================================================================
  ## ============================================================================
  ####### algorithm
  ## ============================================================================
  ## ============================================================================
  ## ============================================================================  
  for(cycle in 1:iters){ # cycle <- 1
    #print(restrained)
    if(!silent){
      count <- count + 1
    }
    
    varos <- lapply(sigma,up.to.vec)
    sigma2 <- as.matrix(unlist(varos)) # current values of vc
    #####################
    ## Expand ZKZ to multitrait
    ## form each T * ZKZ where T is the trait matrix; i.e. sigma[[1]]
    listGs <- list()
    for(i in 1:(nz+nr)){ # for each random effect
      listGs[[i]] <- Matrix::kronecker(sigma[[i]],as.matrix(ZKZ[[i]]))
      # each random effect ZKZ has its own sigma, length of ZKZ and sigma are equal
    }
    #####################
    ### For V and V.inv matrices by adding up all multitrait ZKZ'
    W <- matrix(0,dimos[1]*dimos[2],dimos[1]*dimos[2])
    for(l in 1:length(listGs)){ W <- W + listGs[[l]]};#(W[1:5,1:5]);dim(W)
    
    V <- try(solve(as(W, Class=Vmat.type),sparse=FALSE), silent = TRUE)
    if(class(V) == "try-error"){
      V <- try(solve(as(W + tolparinv * diag(dim(W)[2]), Class=Vmat.type),sparse=FALSE), silent = TRUE)
    }
    W<-NULL
    #####################
    ### PROJECTION MATRIX
    # P = V-  -  V-X [X' V- X]-1 XV- =  WQK
    VX <- V %*% X # V- X
    
    tXVXVX <- try(solve(t(X)%*%VX, t(VX)),silent = TRUE)
    if(class(tXVXVX) == "try-error"){
      tXVXVX <- try(solve((t(X)%*%VX + (tolparinv * diag(dim(t(X)%*%VX)[2]))),t(VX)), silent = TRUE)
    }
    
    P <- V - VX %*% tXVXVX #WQK
    V <- NULL
    # y'Py
    rss <- as.numeric(t(Y) %*% P %*% Y)
    
    sigma2 <- sigma2 *rss/rankX
    coef2[!pos] <- sigma2[!pos] # FALSE are copied
    coef2[pos] <- log(sigma2[pos]) # TRUE are copied
    P <- P * rankX/rss # P * [r(X) / yPy] #WQX <- WQX * rankQK/rss
    rss <- rankX 
    eig <- sort(eigen(P,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankX]
    
    if(any(eig < 0)){
      P <- P + (tolpar - min(eig))*diag(dim(P)[1])
      eig <- eig + tolpar - min(eig)
    }
    
    ldet <- sum(log(eig)) 
    llik <- ldet/2 - rss/2 #.5 [log(det(eigen(P))) - yPy]
    
    if(cycle == 1) llik0 <- llik #keep first log likelihood
    delta.llik <- llik - llik0 # increase of likelihood with respect to initial LL
    llik0 <- llik # update likelihood to current iteration
    
    x <- NULL # a clean x
    var.components <- rep(1,k) # variance components
    ind <- which(pos) # which are TRUE
    if(length(ind)) var.components[ind] <- sigma2[ind] #update
    
    #####################
    ## Multivariate PVi ; where Vi is the derivative ZKZ'
    deriv.list.vc <- list()
    for(v in 1:(nz+nr)){ ## For each random effect obtain all P %*% (T * ZKZ') = P Vi
      deriva <- ZKZ[[v]] # ZKZ
      deriv.list.vc[[v]] <- lapply(posmats.list.vc[[v]],function(x,y){P %*% as(Matrix::kronecker(x,y),Class="sparseMatrix")},y=as.matrix(deriva))# P (T * ZKZ)
    } #tcrossprod(ZETA[[v]]$Z %*% ZETA[[v]]$K, (ZETA[[v]]$Z)) == Vi or dV/ds
    ### convert a 2 level list into a 1 level list
    TT <- do.call(list, unlist(deriv.list.vc, recursive=FALSE))# list of PVi
    #length(TT)
    
    #####################
    ## obtain first derivatives
    # Vi = dV/ds  ..... y' P Vi P y - tr(P Vi)  same than -tr(PVi) - y'PViPy
    x <- sapply(TT,function(x) as.numeric(t(Y) %*% x %*% P %*% Y - sum(diag(x))))
    ## theta(k) * dL/ds  ..... are scalar values
    x <- x * var.components
    
    #####################
    ## second derivatives .... [theta(i) * 1st.deriv(i)] * [theta(j) * 1st.deriv(j)]  * sigma(i) * sigma(j)
    A <- matrix(rep(0, k^2), k, k)
    entries <- expand.grid(1:k,1:k) # indices to be filled
    entries <- entries[!duplicated(t(apply(entries, 1, sort))),] # only for lower triangular
    ## Fisher's Information tr(PVi * PVi) .... A*=Vi=dV/ds .... [Vi Vj'] si sj ; TT is the list of derivatives for all random effects - trait combos
    ff <- function(x) sum(TT[[x[1]]] * t(TT[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
    aa <- apply(entries,1,ff) # matrix of combinations of var.comp 
    A[as.matrix(entries)] <- aa
    A <- copying2(A) # copy lower in upper triangular
    A.svd <- ginv(A) # Inverse of Fishers
    
    #stats <- c(stats, llik, sigma[1:k], x[1:k]) # keep statistics LL, sigma, x * var.components
    ## F- * sigma(k) * dL/ds
    newx <- A.svd %*% x #update 
    #coef is a copy of sigma to be updated
    
    ######################################
    ######################################
    ## constraint to silence zero var comps
    #diagss; traitm
    try1 <- coef2 + taper[cycle] * newx # sigma + f[s*F-*dL/ds] ..... = coef + taper[x]
    if(constraint){
      bad <- diagss2[which(try1[diagss2,] <= 0)]
      #^^^^^ if user wants to set some variance components to zero we will restrain them
      if(!is.null(restrained)){bad <- sort(unique(c(bad,restrained)))}
      
      if(length(bad) > 0){
        #cat("\nRestraining parameters")
        ## find which random effect it is and silence
        badRE <- which(unlist(lapply(diagss, function(x,y){length(which(x%in%y))}, y=bad)) != 0)
        restrain <- list()
        for(b in badRE){# identify diags and off-diags to silence
          bad.diagon <- diagss[[b]][which(diagss[[b]] %in% badRE)] # bad diagonal elements
          # which offdiags correspond to that bad diagonal element
          restrain[[b]] <- sort(unique(which(all2[,1] %in% bad.diagon | all2[,2] %in% bad.diagon)))
          #restrain[[b]] <- sort(c(,offdiagss[[b]])) # varcomp to restrain
        }
        restrain <- unlist(restrain)
        #^^^^^ if user wants to set some variance components to zero we will restrain them
        if(!is.null(restrained)){restrain <- sort(unique(c(restrain,restrained)))}
        
        no.restrain <- setdiff(1:length(try1), restrain)
        A.svd.good <- ginv(A[no.restrain,no.restrain]) # Inverse of Fishers only for good ones
        newx.good <- A.svd.good %*% x[no.restrain] #update 
        newx[no.restrain] <- newx.good
        newx[restrain] <- 0
      }
    }
    ######################################
    ## end of parameter restraining
    ######################################
    
    coef2 <- coef2 + taper[cycle] * newx # sigma + f[s*F-*dL/ds] ..... = coef + taper[x]
    if(constraint){if(length(bad) > 0){coef2[restrain,1] <- 0}} ## make sure you set to zero the bad varcomp
    sigma2[!pos] <- coef2[!pos] # FALSES are replaced
    sigma2[pos] <- exp(coef2[pos]) # TRUES are replaced
    
    ### reaccomodate var.com from vector to list of matrices
    no.var <- lapply(deriv.list.vc,function(x){length(x)})
    #print(no.var)
    sigma3<-cbind(sigma3,sigma2)
    sigmaxxx <- sigma2
    for(r in 1:length(no.var)){
      si <- 1:no.var[[r]]
      newmat <- matrix(NA,ts,ts)
      sq <- upper.tri(newmat)
      diag(sq) <- TRUE
      babas2 <- which(sq,arr.ind = TRUE)
      babas2 <- babas2[ order(babas2[,1], babas2[,2]), ]
      newmat[babas2] <- sigma2[si,1]
      sigma[[r]] <- copying(newmat)
      sigma2 <- matrix(sigma2[-si,])
    }
    sigma <- lapply(sigma, function(x){colnames(x) <- rownames(x) <- namos; return(x)})
    # sigma 2 must finish empty
    # sigma is again filled as the original sigma but with new estimates
    
    var.comp.ret[[cycle]] <- lapply(sigma, function(x,y,z){(x*y)/z},y=base.var,z=sc.var)
    
    if(cycle > 1 & (delta.llik) < tolpar) {#tolpar*10
      if(!silent){
        setTxtProgressBar(pb, ((tot-1)/tot))### keep filling the progress bar
      }
      convergence <- TRUE
      break
    }
    #cat(paste("LL=",llik,"iter=",cycle))
    llstore[cycle] <- llik
    ############ draw
    if(draw){
      layout(matrix(1:2,2,1))
      plot(llstore,type="l",lwd=2,bty="n",col="cadetblue",las=2,xaxt="n",main="LogLikelihood",ylab="LL", xlab="Iteration")
      axis(1,at=1:iters,labels = 1:iters)
      legend("topleft",legend=round(llik,3),bty="n",cex=0.8)
      plot(sigma3[1,],type="l",lwd=2,bty="n",ylim=c(min(sigma3),max(sigma3)),col="white",las=2,xaxt="n",main="Var-Covar components",ylab="scaled var.comp",xlab="Iteration")
      for(u in 1:dim(sigma3)[1]){
        lines(sigma3[u,],col="red",lty=u)
      }
      axis(1,at=1:iters,labels = 1:iters)
    }
    ########## end of draw
    
    ## end of one cycle
    if(!silent){
      setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    
  }
  #############################
  ### END OF CYCLES
  #############################
  #print(var.comp.ret)
  sigma.scaled <- sigma
  good <- which(llstore==max(llstore))
  theta <- var.comp.ret[[good]]
  #theta <- var.comp.ret[[length(var.comp.ret)]]
  ############################################
  ############################################
  ############################################
  ############################################
  ############################################
  
  sigma <- theta 
  if(!is.null(forced)){sigma <- forced}
  #print(sigma)
  #####################
  ## recalculate Vinv with original scale sigma
  listGs <- list()
  for(i in 1:(nz+nr)){listGs[[i]] <- Matrix::kronecker(as.matrix(ZKZ[[i]]),sigma[[i]])}
  W <- matrix(0,dimos[1]*dimos[2],dimos[1]*dimos[2])
  for(l in 1:length(listGs)){ W <- W + listGs[[l]]}#;(W[1:5,1:5]);dim(W)
  V <- try(solve(as(W, Class=Vmat.type),sparse=TRUE), silent = TRUE)
  if(class(V) == "try-error"){
    V <- try(solve(as(W + tolparinv * diag(ncol(W)), Class=Vmat.type),sparse=FALSE), silent = TRUE)
  }
  W <- NULL
  ######################
  ## AIC & BIC
  AIC = as.vector((-2 * llik ) + ( 2 * dim(X)[2]))
  BIC = as.vector((-2 * llik ) + ( log(length(Y)) * dim(X)[2]))
  
  ######################
  ## variance components
  #out1 <- sigma 
  #sigma <- lapply(sigma,function(x){round(x,5)})
  ######################
  ## beta
  
  ## real variance for fixed effects
  V2 <-V
  indexes <- numeric()
  
  # print(ncol(V2))
  # print(ts)
  # print(indexes)
  for(o in 1:ts){
    indexes <- c(indexes,seq(o,ncol(V2),ts))
  }
  V2 <- V2[indexes,indexes]
  XVX2 <- crossprod(X, V2 %*% X)
  XVXi2 <- try(solve(XVX2), silent = TRUE) # variance of fixed effects
  if(class(XVXi2) == "try-error"){
    XVXi2 <- try(solve((XVX2 + tolparinv * diag(dim(XVX)[2]))), silent = TRUE)
  }
  
  XVX <- crossprod(X, V %*% X)
  XVXi <- try(solve(XVX), silent = TRUE) # variance of fixed effects
  if(class(XVXi) == "try-error"){
    XVXi <- try(solve((XVX + tolparinv * diag(dim(XVX)[2]))), silent = TRUE)
  }
  P <- V - V %*% X %*% solve(XVX, crossprod(X, V))
  beta <- (XVXi %*% crossprod((X), V %*% Y.red.noscale.lin)) # (XVX)-XV-Y .... Y.or3 %*% X %*% solve(crossprod(X))#
  beta <- as.matrix(beta)
  
  ######################
  ## residuals and fitted
  XB <- X%*%beta
  fitted.y.good <- matrix(XB, nrow = nrow(Y.red.noscale.ordim), byrow = FALSE); #in a dataframe xb
  ee <- Y.red.noscale.ordim - fitted.y.good
  ee.lin.int <- matrix(t(Y.red.noscale.ordim) - (t(fitted.y.good)), ncol = 1, byrow = FALSE) # residuals transposed and mixed
  Vi.ee <- V %*%  ee.lin.int # V-(y-Xb)' ... nxn %*% linearized(txn) intercaled
  ######################
  ## Multivariate K == T * K
  varvecG <- list()
  for (k in 1:nz) {
    K <-  ZETA[[k]]$K 
    varvecG[[k]] <- Matrix::kronecker(as.matrix(K), as(sigma[[k]],Class="denseMatrix"))
  }
  ######################
  ######################
  ## BLUPs, Var and PEV blups
  ulist <- list()
  Zulist <- list()
  Zulist.ordim <- list()
  Var.u <- list()
  PEV.u <- list()
  Zforvec.list <- list()
  for (a in 1:nz) { # u = GZ'V- (y- XB)  # (30x30) (30x100) (100x100) (100x1)
    
    lev.re <- dim(ZETA[[a]]$Z)[2] # levels of the random effect
    Zforvec <- as(Matrix::kronecker(t(as.matrix(ZETA[[a]]$Z)),diag(ts)), Class="sparseMatrix") # Z' for GZ'
    Zforvec.list[[a]] <-  Zforvec
    Zforvec2 <- as(Matrix::kronecker(as.matrix(ZETA[[a]]$Z),diag(ts)), Class="sparseMatrix") ## Z for Zu
    Zforvec3 <- as(Matrix::kronecker(as.matrix(ZETA.or[[a]]$Z),diag(ts)), Class="sparseMatrix") ## Z.or for Zu.or
    
    ZKforvec <- varvecG[[a]] %*% Zforvec # GZ'
    provi <- ZKforvec %*% Vi.ee # u.hat = GZ'V-(y-Xb) 
    
    ulist[[a]] <- matrix(provi, nrow = lev.re, byrow = TRUE); colnames(ulist[[a]]) <- colnames(Y.red.noscale.ordim) # u
    Var.u[[a]] <- ZKforvec %*% tcrossprod(P, ZKforvec) # var.u.hat = ZGPZ'G ... sigma^4 ZKP ZK
    PEV.u[[a]] <- varvecG[[a]] - Var.u[[a]] #PEV.u.hat = G - ZGPGZ'
    Zulist[[a]] <- matrix(Zforvec2 %*% provi, nrow = nrow(Y.red.noscale.ordim), byrow = TRUE); colnames(Zulist[[a]]) <- colnames(Y.red.noscale.ordim) # Zu reduced
    Zulist.ordim[[a]] <- matrix(Zforvec3 %*% provi, nrow = nrow(Y.or), byrow = TRUE); colnames(Zulist.ordim[[a]]) <- colnames(Y.or) # Zu fitted
    
    Var.u.bytraits <- list()
    PEV.u.bytraits <- list()
    for(h in 1:ts){ # Var.u and PEV.u are mixed for all traits, we need to split them by traits
      name.trait <- traitnames[h]
      rowcol.totake <- seq(h,ncol(Var.u[[a]]),ts)
      Var.u.bytraits[[name.trait]] <- Var.u[[a]][rowcol.totake,rowcol.totake]
      PEV.u.bytraits[[name.trait]] <- PEV.u[[a]][rowcol.totake,rowcol.totake]
    }
    Var.u[[a]] <- Var.u.bytraits
    PEV.u[[a]] <- PEV.u.bytraits
    
    if(EIGEND & a == 1){
      ulist[[a]] <- tUsi %*% matrix(provi, nrow = lev.re, byrow = TRUE); colnames(ulist[[a]]) <- colnames(Y.red.noscale.ordim) # u
      Var.u[[a]] <- lapply(Var.u[[a]], function(x){tUsi %*% tcrossprod(x, tUsi)}) 
      PEV.u[[a]] <- lapply(PEV.u[[a]], function(x){tUsi %*% tcrossprod(x, tUsi)})  # standard errors (SE) for each individual
      # for fitted values
      Zulist[[a]] <- tUsi %*%  matrix(Zforvec2 %*% provi, nrow = nrow(Y.red.noscale.ordim), byrow = TRUE); colnames(Zulist[[a]]) <- colnames(Y.red.noscale.ordim) # Zu reduced
      Zulist.ordim[[a]] <- tUsi %*% matrix(Zforvec3 %*% provi, nrow = nrow(Y.or), byrow = TRUE); colnames(Zulist.ordim[[a]]) <- colnames(Y.or) # Zu fitted
    }
    
    indnames <- colnames(ZETA.or[[a]]$Z)
    if(!is.null(indnames)){
      rownames(ulist[[a]]) <- indnames
    }
  };#(sigma2.u^2) *  ( crossprod(Z%*%K, P)  %*%  (Z%*%K)   )
  names(ulist) <- varosssZ
  names(Zulist.ordim) <- varosssZ
  names(Var.u) <- varosssZ
  names(PEV.u) <- varosssZ
  
  ######################
  ## conditional residuals
  Zu <- 0
  Zu.ordim <- 0
  for(h in 1:length(Zulist)){
    Zu <- Zu + Zulist[[h]]
    Zu.ordim <- Zu.ordim + Zulist.ordim[[h]]
  }
  ee.cond <- Y.red.noscale.ordim - fitted.y.good - Zu # Y - XB - ZU
  ######################
  ## Fisher inverse
  sigma.cov <- ((A.svd) * 2) 
  FI <- (A)/2
  #### convert FI using pos
  FI.c <- matrix(0,dim(FI)[1],dim(FI)[2])
  FI.c <- FI / tcrossprod((sigmaxxx-1)*pos+1)
  sigma.cova <- try(ginv(FI.c),silent=TRUE)
  
  ################################
  ################################
  # PARAMETERS USING ORIGINAL DATA
  ################################
  dado <- lapply(ZETA, function(x){dim(x$Z)})
  ######################
  ## Fitted values and
  ## conditional residuals
  take <- qr$pivot[1:qr$rank]
  rankdeficient <- ncol(X.or) - length(take)
  if(rankdeficient > 0){
    warning(paste("fixed-effect model matrix is rank deficient so dropping", rankdeficient*ts, "columns / coefficients"), call. = FALSE)
    if(ncol(X.or)>1){ # if other than intercept
      X.or <- as.matrix(X.or[,take]) # added to make sure is finished when no full rank
    }
  }
  ncolxor <- ncol(X.or)
  # print(length(take))
  # print(take)
  # print(dim(X.or))
  #X.or <- as.matrix(X.or)
  
  X.or <- do.call("adiag1", rep(list(X.or), ts));
  
  XB.ordim <- matrix(X.or%*%beta, nrow = nrow(Y.or), byrow = FALSE); #in a dataframe xb
  Y.fitted <- XB.ordim + Zu.ordim
  res.ordim <- Y.fitted - Y.or
  
  layout(matrix(1,1,1))
  if(!silent){
    setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  }
  if(is.null(forced)){fofo <- FALSE}else{fofo<- TRUE}
  
  #beta <- matrix(beta,ncol=length(traitnames),nrow=length(fixnames), byrow=FALSE)
  beta <- matrix(beta,ncol=length(traitnames),nrow=ncolxor, byrow=FALSE)
  
  #print(beta)
  
  if(rankdeficient > 0){
    rownames(beta) <- fixnames[take]
  }else{  rownames(beta) <- fixnames}

  colnames(beta) <- traitnames
  
  
  #monitor
  var.mon <- lapply(var.comp.ret,function(x){
    xq <- lapply(x,up.to.vec)
    xq2 <- as.matrix(unlist(xq))
    return(xq2)
  })
  var.mon2 <- do.call(cbind,var.mon)
  monitor <- rbind(llstore[1:ncol(var.mon2)],var.mon2)
  rownames(monitor)[1] <- "llik"
  colnames(monitor) <- paste("iter",1:ncol(monitor))
  ## bring var.comp and fish.inv to non-scale 
  sigma.nonscale <- unlist(lapply(sigma,function(x){up.to.vec(x)}))
  dd <- rep(up.to.vec(base.var),length(sigma))
  ee <- rep(up.to.vec(sc.var),length(sigma))
  fish.inv.nonscale <- ( sigma.cova * (dd%*%t(dd)) ) / ((ee%*%t(ee)))
  
  return(list(var.comp=sigma, V.inv=V, u.hat = ulist , Var.u.hat = Var.u, 
              beta.hat = beta, Var.beta.hat = XVXi2, fish.inv=sigma.cova, 
              fish.inv.nonscale=fish.inv.nonscale,
              PEV.u.hat = PEV.u, residuals=ee, cond.residuals=ee.cond,
              LL=llik, AIC=AIC, BIC=BIC, X=X, Y= Y.red.noscale.ordim,
              #Zforvec= Zforvec.list, varvecG=varvecG,Ylin=Y.red.noscale.lin,
              dimos=dado, sigma.scaled=sigmaxxx, sigma= sigma.nonscale,
              fitted.y=Y.fitted, fitted.u=Zu.ordim, ZETA=ZETA, used.observations=nona,
              method="MNR",random.effs=varosssZ, forced=fofo, convergence=convergence,
              monitor=monitor, restrained=restrained, Zus=Zulist.ordim, 
              res.ordim=res.ordim))
}
