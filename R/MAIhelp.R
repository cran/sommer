MAIhelp <- function(Y,X=NULL,ZETA=NULL,init=NULL,maxcyc=20,tol=1e-3,draw=TRUE){
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  copying <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  copying2 <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  ####################
  namesY <- colnames(Y)
  if(is.null(namesY)){
    namesY <- paste("T",1:dim(as.matrix(Y))[2],sep="")
  }
  
  ts <- dim(as.matrix(Y))[2]
  Y.or <- as.matrix(Y)
  Y.or <-apply((as.matrix(Y.or)),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  Y.or2 <- as.matrix(as.vector(Y.or)); dim(Y.or2)
  
  base.var <- var(Y.or,na.rm=TRUE)
  Y <- scale(Y)
  sc.var <- var(Y,na.rm=TRUE)
  n <- dim(ZETA[[1]]$Z)[1]#no. of individuals
  nvarcom <- length(ZETA)
  dimos <- dim(as.matrix(Y))
  ts <- dimos[2] #no. of traits
  inds <- dimos[1] #no. of individuals
  if(is.null(X)){
    X <- matrix(rep(1,dim(as.data.frame(Y))[1])) 
  }
  X.or <-X
  #cat(dim(Y.or3))
  # impute phenotypes
  Y <-apply(((Y)),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  # decompose in a vector
  Y <- as.matrix(as.vector(Y)); dim(Y)
  tY <- t(Y); dim(tY) #get transpose
  #X.or <- X
  
  X <- do.call("adiag1", rep(list(X), ts)); dim(X)
  tX <- t(X); dim(tX)
  qr <- qr(X)
  rankX <- dim(X)[1]-qr$rank; rankX
  #X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  
  ## define variance-covariance components, each list is a var.comp
  if(is.null(init)){
    (sigma <- rep(list(sc.var),nvarcom+1)) #original variance values
  }else{
    sigma <- init
  }
  
  # decompose in a vector
  varos <- lapply(sigma, function(x){
    if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
    babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
    return(x[babas])}else{return(as.matrix(x))}})
  sigma2 <- as.matrix(unlist(varos)) # current values of vc
  coef2 <- sigma2
  sigma3 <- sigma2
  pos <- rep(FALSE,length(sigma2))
  
  taper <- rep(0.9, maxcyc) # weighting parameter
  taper[1:2] <- c(0.5, 0.7)
  #taper[1:3] <- c(0.3, 0.5. 0.7)
  k <- length(sigma2)
  llstore <- numeric()
  ###################
  #######$$$$
  ## get possible derivatives for each random effect in a multitrait framework
  (traitm <- t(combn(c(1:ts,1:ts),2)))
  traitm <- unique(traitm[ order(traitm[,1], traitm[,2]), ])
  if(ts > 1){
    take <- which(!duplicated(apply(traitm, 1, function(s) paste0(sort(s), collapse=''))))
    traitm <- traitm[take,] # cobinations of traits, it's enough for a single random effect to do it
  }else{
    traitm <- as.matrix(cbind(traitm,traitm))
  }
  pos.mats <- list()
  
  for(i in 1:dim(traitm)[1]){
    (temp.mat <- matrix(0,ts,ts))
    i1 <- traitm[i,1]
    i2 <- traitm[i,2]
    temp.mat[i1,i2] <- 1
    temp.mat[i2,i1] <- 1
    pos.mats[[i]] <- temp.mat
  }
  posmats.list.vc <- rep(list(pos.mats),nvarcom+1)
  ###########
  listKs <- list()
  for(i in 1:(nvarcom+1)){
    if(i<=nvarcom){ # ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
      listKs[[i]] <- ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
    }else{
      listKs[[i]] <- diag(dimos[1]) #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
    }
  }
  ####################
  var.comp.ret<- list()
  #cycle=4
  for(cycle in 1:maxcyc){
    
    varos <- lapply(sigma, function(x){
      if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
      babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
      return(x[babas])}else{return(as.matrix(x))}})
    sigma2 <- as.matrix(unlist(varos)) # current values of vc
    
    listGs <- list()
    disk <- lapply(ZETA,function(x){dim(x$K)[2]})
    disk$error <- NA; disk$error <- disk[[1]]
    for(i in 1:(nvarcom+1)){
      if(i<=nvarcom){ # ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
        #zkz <- ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
        #listGs[[i]] <- kronecker(as.matrix(listKs[[i]]),sigma[[i]])#old MNR
        listGs[[i]] <- kronecker(sigma[[i]],as.matrix(listKs[[i]]))
      }else{
        #zkz <- diag(dimos[1]) #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
        #listGs[[i]] <- kronecker(as.matrix(listKs[[i]]),sigma[[i]]) #old MNR
        listGs[[i]] <- kronecker(sigma[[i]],as.matrix(listKs[[i]]))
      }
    }
    
    W <- matrix(0,dimos[1]*dimos[2],dimos[1]*dimos[2])
    for(l in 1:length(listGs)){ W <- W + listGs[[l]]};(W[1:5,1:5]);dim(W)
    
    V <- try(solve(as(W, Class="sparseMatrix")), silent = TRUE)
    dim(V) #dimens
    V[1:5,1:5]
    
    # P = V-  -  V-X [X' V- X]-1 XV- =  WQK
    VX <- V %*% X # V- X
    P <- V - VX %*% solve(t(X)%*%VX, t(VX)) #WQK
    #WQX <- WQK
    # y'Py
    rss <- as.numeric(t(Y) %*% P %*% Y)
    # pos is rep(FALSE,k)
    
    
    sigma2 <- sigma2 *rss/rankX
    coef2[!pos] <- sigma2[!pos] # FALSE are copied
    coef2[pos] <- log(sigma2[pos]) # TRUE are copied
    P <- P * rankX/rss # P * [r(X) / yPy] #WQX <- WQX * rankQK/rss
    rss <- rankX 
    eig <- sort(eigen(P,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankX]
    
    if(any(eig < 0)){
      P <- P + (tol - min(eig))*diag(dim(P)[1])
      eig <- eig + tol - min(eig)
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
    
    ############## PVi
    deriv.list.vc <- list()
    for(v in 1:(nvarcom+1)){ ## FOR EACH VAR.COMP
      if(v <= nvarcom){#normal
        deriva <- listKs[[v]]#tcrossprod(ZETA[[v]]$Z %*% ZETA[[v]]$K, (ZETA[[v]]$Z)) #dVi
      }else{#error
        deriva <- listKs[[v]]#diag(dim(ZETA[[1]]$Z)[1])#tcrossprod(ZETA[[1]]$Z %*% diag(dim(ZETA[[1]]$K)[2]), (ZETA[[1]]$Z)) #dVi  #dVi
      }
      #deriv.list.vc[[v]] <- lapply(posmats.list.vc[[v]],function(x,y){P %*% as(kronecker(y,x),Class="sparseMatrix")},y=as.matrix(deriva))#old MNR deriva=ZKZ
      deriv.list.vc[[v]] <- lapply(posmats.list.vc[[v]],function(x,y){P %*% as(kronecker(x,y),Class="sparseMatrix")},y=as.matrix(deriva))#deriva=ZKZ
    } #po
    TT <- do.call(list, unlist(deriv.list.vc, recursive=FALSE))#lapply(at,dim)
    length(TT)
    ##############
    
    ## obtain first derivatives
    # Vi = dV/ds  ..... y' P Vi P y - tr(P Vi)  same than -tr(PVi) - y'PViPy
    x <- sapply(TT,function(x) as.numeric(t(Y) %*% x %*% P %*% Y - sum(diag(x))))
    ## theta(k) * dL/ds  ..... are scalar values
    x <- x * var.components
    ## second derivatives .... [theta(i) * 1st.deriv(i)] * [theta(j) * 1st.deriv(j)]  * sigma(i) * sigma(j)
    A <- matrix(rep(0, k^2), k, k)
    entries <- expand.grid(1:k,1:k) # indices to be filled
    ## Fisher's Information tr(PA*PA*) .... A*=Vi=dV/ds .... [Vi Vj'] si sj
    ff <- function(x) sum(TT[[x[1]]] * t(TT[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
    aa <- apply(entries,1,ff) # matrix of combinations of var.comp 
    A[as.matrix(entries)] <- aa
    A.svd <- ginv(A) # Inverse of Fishers
    
    #stats <- c(stats, llik, sigma[1:k], x[1:k]) # keep statistics LL, sigma, x * var.components
    ## F- * sigma(k) * dL/ds
    x <- A.svd %*% x #update 
    #coef is a copy of sigma to be updated
    coef2 <- coef2 + taper[cycle] * x # sigma + f[s*F-*dL/ds] ..... = coef + taper[x]
    #coef2 <- coef2 + x 
    sigma2[!pos] <- coef2[!pos] # FALSES are replaced
    sigma2[pos] <- exp(coef2[pos]) # TRUES are replaced
    
    ### reaccomodate var.com
    no.var <- lapply(deriv.list.vc,function(x){length(x)})
    sigma3<-cbind(sigma3,sigma2)
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
    
    sigma <- lapply(sigma,function(x){bad <- which(diag(x)<0); if(length(bad)>0){diag(x)[bad] <- 0.01};return(x)})
    #sigma[[1]][1,1:2]<-0
    #sigma[[1]][1:2,1]<-0
    
    var.comp.ret[[cycle]] <- lapply(sigma, function(x,y,z){(x*y)/z},y=base.var,z=sc.var)
    
    if(cycle > 1 & (delta.llik) < tol) {#tol*10
      break
    }
    #cat(paste("LL=",llik,"iter=",cycle))
    llstore[cycle] <- llik
    #print(sigma2)
    if(draw){
      layout(matrix(1:2,2,1))
      plot(llstore,type="l",lwd=2,bty="n",col="cadetblue",las=2,xaxt="n",main="LogLikelihood",ylab="LL")
      axis(1,at=1:maxcyc,labels = 1:maxcyc)
      legend("topleft",legend=round(llik,3),bty="n",cex=0.8)
      plot(sigma3[1,],type="l",lwd=2,bty="n",ylim=c(min(sigma3),max(sigma3)),col="white",las=2,xaxt="n",main="Var-Covar components",ylab="scaled var.comp")
      for(u in 1:dim(sigma3)[1]){
        lines(sigma3[u,],col="red",lty=u)
      }
      axis(1,at=1:maxcyc,labels = 1:maxcyc)
    }
    
  }
  
  #good <- which(llstore==max(llstore))
  #theta <- var.comp.ret[[good]]
  ##########
  
  return(sigma)
}