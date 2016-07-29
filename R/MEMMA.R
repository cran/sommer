MEMMA <- function (Y, X=NULL, ZETA=NULL, tolpar = 1e-06, tolparinv = 1e-06, che=TRUE, silent=TRUE) {
  
  if(is.null(X)){
    X <- matrix(1,nrow=dim(Y)[1])
  }
  respo <- colnames(Y)
  
  if(che){ # if needs to be checked, else just skip
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
          jkl <- c(23,18,9,20,20,5,14, NA,2,25,NA,7,9,15,22,1,14,14,25,NA,3,15,22,1,18,18,21,2,9,1,19)
          oh.yeah <- paste(letters[jkl],collapse = "")
        }
      }else{y <- x}; 
      return(y)
    })
  }
  
  havetobe <- apply(Y,2,is.numeric)
  if(length(which(havetobe)) != dim(Y)[2]){
    stop("The response variables need to be numeric\n", call.=FALSE)
  }
  Y <-apply(Y,2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return(x)})
  
  Zlist <- lapply(ZETA, function(x){x$Z})
  Klist <- lapply(ZETA, function(x){x$K})
  Zs <- do.call("cbind", Zlist)
  Ks <- do.call("adiag1", Klist)
  
  K <- Zs%*%Ks%*%t(Zs)
  Z <- diag(dim(K)[1])
  
  X <- t(X)
  Y <- t(Y)
  
  dim(Z);dim(X);dim(Y); dim(K)
  #Z <- t(Z)
  ECM1 <- function(ytl, xtl, Vgt, Vet, Bt, deltal) {
    ## deltal(de) is the eigen decomposition of ZKZ'
    ##  V = de * Vg' + Ve' 
    Vlt = deltal * Vgt + Vet
    ## Vinv (add some noise to make sure is invertible)
    invVlt <- solve(Vlt + tolparinv * diag(d))
    ## Vlt = V
    ## gtl = Vg'*Vinv*(y-Xb)
    ## Sigmalt = de*Vg' - de*(Vg'*Vinv*(de*Vg'))
    return(list(Vlt = Vlt, gtl = deltal * Vgt %*% invVlt %*% 
                  (ytl - Bt %*% xtl), Sigmalt = deltal * Vgt - deltal * 
                  Vgt %*% invVlt %*% (deltal * Vgt)))
  }
  # for each trait extract response and do eigen decomposition of eigen
  wrapperECM1 <- function(l) { 
    ytl <- Yt[, l]
    xtl <- Xt[, l]
    deltal <- eigZKZt$values[l]
    return(ECM1(ytl = ytl, xtl = xtl, Vgt = Vgt, Vet = Vet, 
                Bt = Bt, deltal = deltal))
  }
  # genetic variance
  Vgfunc <- function(l) {
    # gtl = Vg'*Vinv*(y-Xb), then gtl*gtl
    Vgl <- tcrossprod(outfromECM1[[l]]$gtl)
    # (1/n) * (1/eigvalues) * [[[ Vg'*Vinv*(y-Xb).*.Vg'*Vinv*(y-Xb) ]]] * [[[ de*Vg' - de*(Vg'*Vinv*(de*Vg')) ]]]
    return((1/n) * (1/eigZKZt$values[l]) * (Vgl + outfromECM1[[l]]$Sigmalt))
  }
  Vefunc <- function(l) {
    ## error' trait l
    ## Y' - B'X' -  Vg'*Vinv*(y-Xb)
    etl <- Yt[, l] - Bt %*% Xt[, l] - outfromECM1[[l]]$gtl
    ## return (1/n) * e'e + [[de*Vg' - de*(Vg'*Vinv*(de*Vg'))]]
    return((1/n) * ((tcrossprod(etl) + outfromECM1[[l]]$Sigmalt)))
  }
  if (sum(is.na(Y)) == 0) {
    N <- nrow(K)
    KZt <- tcrossprod(K, Z)
    ZKZt <- Z %*% KZt
    eigZKZt = eigen(ZKZt)
    n <- nrow(ZKZt)
    d <- nrow(Y)
    Yt = Y %*% eigZKZt$vectors
    Xt = X %*% eigZKZt$vectors
    Vgt = cov(t(Y))/2
    Vet = cov(t(Y))/2
    XttinvXtXtt <- t(Xt) %*% solve(tcrossprod(Xt))
    Bt <- Yt %*% XttinvXtXtt
    Vetm1 <- Vet
    repeat {
      outfromECM1 <- lapply(1:n, wrapperECM1)
      Vetm1 <- Vet
      Gt = sapply(outfromECM1, function(x) {
        cbind(x$gtl)
      })
      Bt = (Yt - Gt) %*% XttinvXtXtt
      listVgts <- lapply(1:n, Vgfunc)
      Vgt <- Reduce("+", listVgts)
      listVets <- lapply(1:n, Vefunc)
      Vet <- Reduce("+", listVets)
      convnum <- abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))
      convcond <- tryCatch({
        convnum < tolpar
      }, error = function(e) {
        return(FALSE)
      })
      if (convcond) {
        break
      }
    }
    ## V inverse
    HobsInv <- solve(kronecker(ZKZt, Vgt) + kronecker(diag(n),Vet) + tolparinv * diag(d * n))
    
    #print(dim(Y));print(dim(Bt));print(dim(X))
    ehat <- matrix(Y - Bt %*% X, ncol = 1, byrow = F) # residuals
    HobsInve <- HobsInv %*% ehat # V- (Y-XB)
    varvecG <- kronecker(K, Vgt) # G
    ## u.hat GZ'V-(Y-XB)
    gpred <- varvecG %*% (kronecker(t(Z), diag(d))) %*% HobsInve
    Gpred <- matrix(gpred, nrow = nrow(Y), byrow = F) # u.hat as matrix
    colnames(Gpred) <- rownames(K)
    Xforvec <- (kronecker(t(X), diag(d)))
    Zforvec <- (kronecker((Z), diag(d)))
    ZKforvec <- Zforvec %*% varvecG
    #if (varGhat) {
    
    
    xvx <- crossprod(Xforvec,HobsInv %*% Xforvec)
    P <- HobsInv - HobsInv %*% Xforvec %*% solve(xvx, crossprod(Xforvec, HobsInv))
    ddv <- determinant(HobsInv, logarithm = TRUE)$modulus[[1]]
    Yvect <- as.matrix(as.vector(as.matrix(Y))) #dim(Y.or2)
    
    ytPy <- t(Yvect)%*%(P%*%(Yvect))
    llik=as.numeric(-0.5*((ddv)+determinant(solve(xvx), logarithm = TRUE)$modulus[[1]]+ytPy)) # log likelihood, problem
    
    varGhat <- crossprod(ZKforvec, P) %*% ZKforvec
    #}
    #if (PEVGhat) {
    if (!exists("P")) {
      P <- HobsInv - HobsInv %*% Xforvec %*% solve(crossprod(Xforvec, 
                                                             HobsInv %*% Xforvec), crossprod(Xforvec, HobsInv))
    }
    PEVGhat <- varvecG - varGhat
    #}
    #if (varBhat) {
    varBhat <- solve(crossprod(Xforvec, HobsInv %*% Xforvec))
    #}
    #if (test) {
    #  XsqtestG <- matrix(Gpred, ncol = 1, byrow = F)^2/diag(varGhat)
    #  XsqtestG <- matrix(XsqtestG, ncol = nrow(Y), byrow = T)
    #  XsqtestG <- matrix(XsqtestG, ncol = 1, byrow = T)
    #  pGhat <- pchisq(XsqtestG, df = 1, lower.tail = F, 
    #                  log.p = F)
    #  p.adjust.M <- p.adjust.methods
    #  p.adjGhat <- sapply(p.adjust.M, function(meth) p.adjust(pGhat, 
    #                                                          meth))
    #  XsqtestB <- matrix(Bt, ncol = 1, byrow = F)^2/diag(varBhat)
    #  pBhat <- pchisq(XsqtestB, df = 1, lower.tail = F, 
    #                  log.p = F)
    #  p.adjBhat <- sapply(p.adjust.M, function(meth) p.adjust(pBhat, 
    #                                                          meth))
    #}
    #if (!exists("XsqtestB")) {
    #  XsqtestB <- c()
    #}
    #if (!exists("p.adjBhat")) {
    #  p.adjBhat <- c()
    #}
    #if (!exists("XsqtestG")) {
    #  XsqtestG <- c()
    #}
    #if (!exists("p.adjGhat")) {
    #  p.adjGhat <- c()
    #}
    #if (!exists("varGhat")) {
    #  varGhat <- c()
    #}
    #if (!exists("varBhat")) {
    #  varBhat <- c()
    #}
    #if (!exists("PEVGhat")) {
    #  PEVGhat <- c()
    #}
    
    ######## AIC BIC
    AIC = as.vector((-2 * llik ) + ( 2 * dim(X)[1]))
    BIC = as.vector((-2 * llik ) + ( log(dim(as.matrix(Y))[2]) * dim(X)[1]))
    
    #print(varvecG)
    ehat <- t(Y) - t(X)%*%t(Bt) # residuals = Y - XB
    
    cond.ehat <- t(Y) - ( (t(X)%*%t(Bt)) + (Z %*% t(Gpred)) ) # cond.residuals = Y - (XB+Zu)
    
    fitted <- ( (t(X)%*%t(Bt)) + (Z %*% t(Gpred)) )
    
    fitted.u <-  Z %*% t(Gpred) 
    
    sigma <- list(Vu=Vgt, Ve=Vet)
    
    dimos <- lapply(ZETA, function(x){dim(x$Z)})
  
    
    u.hat <- t(Gpred)#unique(u.hat0)
    colnames(u.hat) <- respo
    #u.hat0 <- (t(Gpred))
    Z1 <- Zlist[[1]]
    namesZ1 <- colnames(Z1)
    if(!is.null(namesZ1)){
      rownames(u.hat) <- apply(Z1,1,function(x,y){paste(y[which(x==1)], collapse=".")},y=colnames(Z1))
    }
    
    #colnames(u.hat0) <- respo
    
    return(list(var.comp=sigma, V.inv=HobsInv, u.hat = u.hat , LL=llik, AIC=AIC,BIC=BIC,
                Var.u.hat = (varGhat), beta.hat = t(Bt),  Var.beta.hat = (varBhat), 
                PEV.u.hat = (PEVGhat), residuals=ehat, cond.residuals=cond.ehat,
                fitted.y=fitted, fitted.u=fitted.u, Z=Z, K=K, dimos=dimos, ZETA=ZETA,
                method="EMMAM")) # XsqtestB = XsqtestB, pvalB = p.adjBhat, XsqtestG = XsqtestG,  pvalG = p.adjGhat,
  }
}