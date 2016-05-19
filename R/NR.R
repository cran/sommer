# Newton-Raphson
NR <- function(y, X=NULL, ZETA=NULL, R=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=50, constraint=TRUE, init=NULL, sherman=FALSE, che=TRUE, MTG2=FALSE, Fishers=FALSE, gss=TRUE, forced=NULL){
  y.or <- y
  x.or <- X
  ### make full function
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
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
  
  #df.ord <- 1:length(ZETA)
  
  ###########################
  ## y is a vector for the response variable
  ## X is an incidence matrix for fixed effects
  ## Z is a list of lists for each random effect
  # the list of Z can or cannot include the covariance matrix for such random effect
  # if provided must be provided as Z=list(list(Z=Z1,K=K1),list(Z=Z2,K=K2), etc) 
  ############################
  if(che){ # if coming from mmer don't check
    if(is.list(ZETA)){
      if(is.list(ZETA[[1]])){ # if was provided as a two level list
        ZETA=ZETA
      }else{ # if was provided as a one level list
        ZETA=list(ZETA)
      }
    }else{
      #stop;
      cat("\nThe random effects need to be provided in a list format, please see examples")
    }
  }
  ###########################
  # if X matrix is not present
  if(is.null(X) & is.null(ZETA)){ # nothing in the model
    tn = length(y); xm <- matrix(1,tn,1)
    yv <- y#scale(y)
    res <- lm(yv~xm-1) # intercept model
  }else{
    if(is.null(X) & !is.null(ZETA)){ # only random effects present
      tn = length(y); xm <- matrix(1,tn,1)
    }
    if(!is.null(X) & !is.null(ZETA)){ # both present, extract xm from X, check double list
      if(is.list(X)){
        if(is.list(X[[1]])){
          xm=X[[1]][[1]]
        }else{
          xm=X[[1]]
        }
      }else{
        xm=as.matrix(X) 
      }
    }
    ############################################
    ## if K matrices are not present in ZETA
    # add an identity matrix to all effects in ETA that did not provide a var-cov matrix
    if(is.null(R)){R <- diag(length(y))} # dim(x[[1]])[2]
    
    if(che){ # if needs to be checked, else just skip
      ZETA <- lapply(ZETA, function(x){
        if(length(x) == 1){
          provided <- names(x)
          if(provided == "Z"){
            y <- list(Z=x[[1]],K=diag(dim(x[[1]])[2]))
          }
          if(provided == "K"){
            y <- list(Z=diag(length(y)),K=x[[1]])
          }else{
            stop()
            cat("Names of matrices provided can only be 'Z' or 'K', the names you provided don't match the arguments required")
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
    if(!is.null(init)){init <- init[c(df.ord, (length(init)))]}
    if(!is.null(forced)){forced <- forced[c(df.ord, (length(forced)))]}
    #####################################################
    ## to use later for fitted values
    x.or <- as.matrix(xm)
    zeta.or <- ZETA
    zeta.or  <- lapply(zeta.or , function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    ##
    if(length(ZETA)==1 & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
      misso <- which(is.na(y))
      if(length(misso) >0){
        y[misso] <- median(y, na.rm=TRUE)
      }
    }
    ZETA2 <- ZETA; y2 <- y ; good <- which(!is.na(y)) # make a copy and identify no missing data
    #ZETA <- lapply(ZETA2, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
    if(length(ZETA)==1 & MTG2==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
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
    xm <- as.matrix(xm[good,])
    txm <- t(xm)
    R <- R[good,good]
    
    Zs <- lapply(ZETA, function(x){x[[1]]})
    Gs <- lapply(ZETA, function(x){x[[2]]})
    Zsp <- as(do.call("cbind", Zs),Class="sparseMatrix") # column bind Z=[Z1 Z2 Z3]
    tZsp <- t(Zsp)
    Ksp <- as(do.call("adiag1", Gs),Class="sparseMatrix") # G as diagonal
    
    ##################################
    ## CREATE V LIST WITH ZKZ MATRICES
    ##################################
    V <- list()
    for(h in 1:(length(ZETA)+1)){
      if(h <= length(ZETA)){
        V[[h]] <- tcrossprod(ZETA[[h]][[1]], ZETA[[h]][[1]] %*% (ZETA[[h]][[2]]) )  
      }else{
        V[[h]] <- R
      }
    }
    ##################################
    ## START NEWTON-RAPHSON ALGORITHM
    ##################################
    if(is.null(names(ZETA))){
      varosss <- c(paste("u.",df.ord, sep=""))
    }else{
      varosss <- as.character(names(ZETA))
    }
    
    k <- length(V)
    n <- length(y)
    maxcyc <- iters
    tol <- 1e-4
    delta <- 1
    pos <- rep(FALSE,k)
    SWsolveINDICATOR <- FALSE
    In <- R
    K <- xm
    qr <- qr(K)
    rankQK <- n - qr$rank
    logL2.stored <- numeric()
    wi=0 # counter
    
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
    #k <- length(V)
    A <- matrix(rep(0, k^2), k, k)
    entries <- expand.grid(1:k,1:k)
    x <- rep(0,k)
    sigma <- c(1,rep(0, k-1))
    stats <- rep(0, 0)
    
    # taper is the steps to take 0-1
    taper <- rep(0.9, maxcyc)
    # starting values for the random effects
    if(!is.null(init)) {
      start <- init[1:k]
    }else{
      start <- rep(var(y,na.rm=TRUE),k)
    }
    #########
    sigma <- coef <- start
    record <- as.matrix(sigma)
    ## reparameterise so everything will get into the correct spot after exp
    coef[pos] <- log(sigma[pos])
    coef[!pos] <- sigma[!pos]
    
    ## Set the memory requirements beforehand
    TT <- vector("list", length=k)
    for(ii in 1:k){TT[[ii]] <- matrix(NA,n,n)}
    
    ##################
    ## progress bar
    if(!silent){
      count <- 0
      tot <- 15
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    ##################
    if(is.null(forced)){ # &&&&&&&&&&&&&&&& IF NOT FORCED &&&&&&&&&&&&&&&&&
      for(cycle in 1:maxcyc){
        
        if(!silent){
          count <- count + 1
        }
        wi <- wi+1
        ## Limit how far out we go on the logarithmic scale
        ind <- which(pos)
        if(length(ind)) {
          coef[ind] <- pmin(coef[ind],20)
          coef[ind] <- pmax(coef[ind],-20) ## so on regular scale everything is between exp(-20) and exp(20)
          sigma[ind] <- exp(coef[ind])
        }
        
        #if(verbose) {
        #  cat(cycle, "sigma =",sigma)
        ##cat(sigma)
        #}
        
        if(!SWsolveINDICATOR) {
          Sigma <- 0
          ## can we get rid of this loop? multiplies each var-cov for its variance component
          ## additionally sums up each var-cov V = ZKZ + ZKZ + ... + ZRZ
          for(i in 1:k) Sigma <- Sigma + V[[i]]*sigma[i]
          # solves V to get V.inverse using In (identity matrix) as response
          W <- solve(Sigma,In)
        } else {
          W <- SWsolve2(Z[1:(k-1)],sigma)
        }
        #### K is the X matrix for fixed effects
        #if(is.null(K)){
        #  WQK <- W
        #}else{ # our case, use sherman
        WK <- W %*% K # V.inv X
        WQK <- W - WK %*% solve(t(K)%*%WK, t(WK)) # V.inv - V.inv X [X' V.inv X]-1
        #}
        ## WQK and WQX
        if(REML){ # REML
          WQX <- WQK
        }else{ # ML
          KX <- K
          WX <- W %*% KX        # including the kernel (Oct 12 2011)
          WQX <- W - WX %*% solve(t(KX)%*%WX, t(WX))
        }
        
        # y' [parameter space] y
        rss <- as.numeric(t(y) %*% WQX %*% y)
        
        sigma <- sigma * rss/rankQK
        coef[!pos] <- sigma[!pos]
        coef[pos] <- log(sigma[pos])
        WQK <- WQK * rankQK/rss
        WQX <- WQX * rankQK/rss
        rss <- rankQK ## looks bad but the rss is absorbed into WQK so the rss term comes out of eig below
        
        eig <- sort(eigen(WQK,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankQK]
        if(any(eig < 0)){
          cat("error: Sigma is not positive definite on contrasts: range(eig)=", range(eig), "\n")
          WQK <- WQK + (tol - min(eig))*diag(nobs)
          eig <- eig + tol - min(eig)
        }
        ldet <- sum(log(eig))
        llik <- ldet/2 - rss/2
        if(cycle == 1) llik0 <- llik
        delta.llik <- llik - llik0
        llik0 <- llik
        logL2.stored <- c(logL2.stored,llik0)
        ## From Jean-Luc Jannick to fix excess carriage returns
        #if (verbose){
        #  if (reml) cat(" resid llik =", llik, "\n") else cat(" llik =", llik, "\n")
        #  cat(cycle, "adjusted sigma =", sigma)
        #  if (cycle > 1){
        #    if (reml) cat(" delta.llik =", delta.llik, "\n") else cat(" delta.llik =", delta.llik, "\n")
        #  } else cat("\n")
        #}
        
        ## now the fun starts, derivative and expected fisher info
        ## the 0.5 multiple is ignored, it is in both and they cancel
        
        ##T <- list(NULL)
        x <- NULL
        identity=TRUE
        
        ## derivatives are now D[[i]] = var.components[i]*V[[i]]
        var.components <- rep(1,k)
        ind <- which(pos)
        if(length(ind)) var.components[ind] <- sigma[ind]
        
        ## Slow part - order k n-squared
        if(!SWsolveINDICATOR) {
          if(identity) {
            TT[[k]] <- WQK
            if(k>1) {
              for(ii in (k-1):1) TT[[ii]] <- WQK %*% V[[ii]] ## # V.inv - V.inv X [X' V.inv X]-1 * ZKZi
            }
          } else {
            for(ii in 1:k) TT[[ii]] <- WQK %*% V[[ii]]
          }
        } else {
          if(identity) {
            TT[[k]] <- WQK
            if(k>1) {
              for(ii in (k-1):1) TT[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
            }
          } else {
            for(ii in 1:k) TT[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
          }
        }
        
        # V.inv - V.inv X [X' V.inv X]-1, adjust +- variance components
        x <- sapply(TT,function(x) as.numeric(t(y) %*% x %*% WQX %*% y - sum(diag(x))))
        x <- x * var.components
        
        ## See nested for loops commented out below - evaluating the Expected Fisher Information, A
        ff <- function(x) {sum(TT[[x[1]]] * t(TT[[x[2]]])) * var.components[x[1]] * var.components[x[2]]}
        aa <- apply(entries,1,ff)
        # A is the Fishers information matrix
        A[as.matrix(entries)] <- aa
        
        stats <- c(stats, llik, sigma[1:k], x[1:k])
        #if(verbose>=9) {
        ##cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
        #}
        
        A.svd <- ginv(A)
        x <- A.svd %*% x
        
        if(qr(A)$rank < k){
          if(cycle==1) {
            #if(verbose) {
              cat("Warning: Non identifiable dispersion model\n")
              ##print(round(A,6))
              cat(sigma)
              cat("\n")
            #}
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
        record <- cbind(record, as.matrix(sigma))
        #############
        ## IF USER WANT TO SEE THE LIKELIHOOD PLOT
        if(draw){
          ylim <- max(unlist(record), na.rm=TRUE)
          my.palette <- brewer.pal(7,"Accent")
          layout(matrix(1:2,2,1))
          plot(logL2.stored,type="l", main="logLikelihood", col=my.palette[7],lwd=3, las=2, xaxt="n", ylab="logLikelihood value", xlab="Iterations processed", cex.axis=0.5) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          legend("bottomleft", legend = round(llik0,3), bty="n", cex=0.7)
          plot(record[1,],ylim=c(0,ylim),type="l", las=2, xaxt="n",main="Newton-Raphson algorithm results", col=my.palette[1],lwd=3, ylab="Value of the variance component", xlab="Iterations processed", cex.axis=0.5) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          for(t in 1:(dim(record)[1])){
            lines(record[t,],col=my.palette[t],lwd=3)
          } 
          
          ww <- dim(record)[1]
          lege <- list()
          #lege2 <- list()
          for(k in 1:length(V)){
            if(k == length(V)){
              lege[[k]] <- paste("Var(e):",round(record[k,wi+1],4), sep="")
              #lege2[[k]] <- paste("Var(e):")
            }else{
              lege[[k]] <- paste("Var(",varosss[k],"):",round(record[k,wi+1],4), sep="")
              #lege2[[k]] <- paste("Var(u",k,"):",sep="")
            }
          }
          legend("topleft",bty="n", cex=0.7, col=my.palette, lty=1, lwd=3, legend=unlist(lege))
        }
        ############# END OF LIKELIHOOD PLOT
        
        ## keep filling the progressbar
        if(!silent){
          setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
        }
        ##############
        
        ## check the change in llik is small
        if(cycle > 1 & abs(delta.llik) < tol*10){
          if(!silent){
            setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
          }
          break
        } 
        if(max(abs(x)) < tol){
          if(!silent){
            setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
          }
          break
        }
      }##################################
      ## END OF NEWTON-RAPHSON ALGORITHM
      ##################################
      badd <- which(sigma < 0 )
      goodd <- which(sigma > 0 )
      ###*****************************************************************************************
      if(constraint & length(badd) > 0){ # if a variance component was negative and constraint was activated
        # recalculate variance components removing such variance component
        cat("\nOne or more variance components close to zero. Boundary constraint applied.\n")
        ZETAXXX <- ZETA[-badd]
        ##################################
        ## CREATE V LIST WITH ZKZ MATRICES
        ##################################
        V <- list()
        for(h in 1:(length(ZETAXXX)+1)){
          if(h <= length(ZETAXXX)){
            V[[h]] <- tcrossprod(ZETAXXX[[h]][[1]], ZETAXXX[[h]][[1]] %*% (ZETAXXX[[h]][[2]]) )  
          }else{
            V[[h]] <- R
          }
        }
        ZETAXXX <- NULL
        ##################################
        ## START NEWTON-RAPHSON ALGORITHM
        ##################################
        k <- length(V)
        n <- length(y)
        maxcyc <- iters
        tol <- 1e-4
        delta <- 1
        pos <- rep(FALSE,k)
        SWsolveINDICATOR <- FALSE
        In <- R
        K <- xm
        qr <- qr(K)
        rankQK <- n - qr$rank
        logL2.stored <- numeric()
        wi=0 # counter
        
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
        #k <- length(V)
        A <- matrix(rep(0, k^2), k, k)
        entries <- expand.grid(1:k,1:k)
        x <- rep(0,k)
        sigma <- c(1,rep(0, k-1))
        stats <- rep(0, 0)
        
        # taper is the steps to take 0-1
        taper <- rep(0.9, maxcyc)
        # starting values for the random effects
        if(!is.null(init)) {
          start <- init[1:k]
        }else{
          start <- rep(var(y,na.rm=TRUE),k)
        }
        #########
        sigma <- coef <- start
        record <- as.matrix(sigma)
        ## reparameterise so everything will get into the correct spot after exp
        coef[pos] <- log(sigma[pos])
        coef[!pos] <- sigma[!pos]
        
        ## Set the memory requirements beforehand
        TT <- vector("list", length=k)
        for(ii in 1:k){TT[[ii]] <- matrix(NA,n,n)}
        
        ##################
        ## progress bar
        if(!silent){
          count <- 0
          tot <- iters
          pb <- txtProgressBar(style = 3)
          setTxtProgressBar(pb, 0)
        }
        ##################
        
        for(cycle in 1:maxcyc){
          
          if(!silent){
            count <- count + 1
          }
          wi <- wi+1
          ## Limit how far out we go on the logarithmic scale
          ind <- which(pos)
          if(length(ind)) {
            coef[ind] <- pmin(coef[ind],20)
            coef[ind] <- pmax(coef[ind],-20) ## so on regular scale everything is between exp(-20) and exp(20)
            sigma[ind] <- exp(coef[ind])
          }
          
          #if(verbose) {
          #  cat(cycle, "sigma =",sigma)
          ##cat(sigma)
          #}
          
          if(!SWsolveINDICATOR) {
            Sigma <- 0
            ## can we get rid of this loop? multiplies each var-cov for its variance component
            ## additionally sums up each var-cov V = ZKZ + ZKZ + ... + ZRZ
            for(i in 1:k) Sigma <- Sigma + V[[i]]*sigma[i]
            # solves V to get V.inverse using In (identity matrix) as response
            W <- solve(Sigma,In)
          } else {
            W <- SWsolve2(Z[1:(k-1)],sigma)
          }
          #### K is the X matrix for fixed effects
          #if(is.null(K)){
          #  WQK <- W
          #}else{ # our case, use sherman
          WK <- W %*% K # V.inv X
          WQK <- W - WK %*% solve(t(K)%*%WK, t(WK)) # V.inv - V.inv X [X' V.inv X]-1
          #}
          ## WQK and WQX
          if(REML){ # REML
            WQX <- WQK
          }else{ # ML
            WX <- W %*% KX        # including the kernel (Oct 12 2011)
            WQX <- W - WX %*% solve(t(KX)%*%WX, t(WX))
          }
          
          # y' [parameter space] y
          rss <- as.numeric(t(y) %*% WQX %*% y)
          
          sigma <- sigma * rss/rankQK
          coef[!pos] <- sigma[!pos]
          coef[pos] <- log(sigma[pos])
          WQK <- WQK * rankQK/rss
          WQX <- WQX * rankQK/rss
          rss <- rankQK ## looks bad but the rss is absorbed into WQK so the rss term comes out of eig below
          
          eig <- sort(eigen(WQK,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankQK]
          if(any(eig < 0)){
            cat("error: Sigma is not positive definite on contrasts: range(eig)=", range(eig), "\n")
            WQK <- WQK + (tol - min(eig))*diag(nobs)
            eig <- eig + tol - min(eig)
          }
          ldet <- sum(log(eig))
          llik <- ldet/2 - rss/2
          if(cycle == 1) llik0 <- llik
          delta.llik <- llik - llik0
          llik0 <- llik
          logL2.stored <- c(logL2.stored,llik0)
          ## From Jean-Luc Jannick to fix excess carriage returns
          
          ## now the fun starts, derivative and expected fisher info
          ## the 0.5 multiple is ignored, it is in both and they cancel
          
          ##T <- list(NULL)
          x <- NULL
          identity=TRUE
          
          ## derivatives are now D[[i]] = var.components[i]*V[[i]]
          var.components <- rep(1,k)
          ind <- which(pos)
          if(length(ind)) var.components[ind] <- sigma[ind]
          
          ## Slow part - order k n-squared
          if(!SWsolveINDICATOR) {
            if(identity) {
              TT[[k]] <- WQK
              if(k>1) {
                for(ii in (k-1):1) TT[[ii]] <- WQK %*% V[[ii]] ## # V.inv - V.inv X [X' V.inv X]-1 * ZKZi
              }
            } else {
              for(ii in 1:k) TT[[ii]] <- WQK %*% V[[ii]]
            }
          } else {
            if(identity) {
              TT[[k]] <- WQK
              if(k>1) {
                for(ii in (k-1):1) TT[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
              }
            } else {
              for(ii in 1:k) TT[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
            }
          }
          
          # V.inv - V.inv X [X' V.inv X]-1, adjust +- variance components
          x <- sapply(TT,function(x) as.numeric(t(y) %*% x %*% WQX %*% y - sum(diag(x))))
          x <- x * var.components
          
          ## See nested for loops commented out below - evaluating the Expected Fisher Information, A
          ff <- function(x) {sum(TT[[x[1]]] * t(TT[[x[2]]])) * var.components[x[1]] * var.components[x[2]]}
          aa <- apply(entries,1,ff)
          # A is the Fishers information matrix
          A[as.matrix(entries)] <- aa
          
          stats <- c(stats, llik, sigma[1:k], x[1:k])
          #if(verbose>=9) {
          ##cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
          #}
          
          A.svd <- ginv(A)
          x <- A.svd %*% x
          
          if(qr(A)$rank < k){
            if(cycle==1) {
              #if(verbose) {
                cat("Warning: Non identifiable dispersion model\n")
                ##print(round(A,6))
                cat(sigma)
                cat("\n")
              #}
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
          record <- cbind(record, as.matrix(sigma))
          #############
          ## IF USER WANT TO SEE THE LIKELIHOOD PLOT
          if(draw){
            ylim <- max(unlist(record), na.rm=TRUE)
            my.palette <- brewer.pal(7,"Accent")
            layout(matrix(1:2,2,1))
            plot(logL2.stored,type="l", main="logLikelihood", col=my.palette[7],lwd=3, las=2, xaxt="n", ylab="logLikelihood value", xlab="Iterations processed", cex.axis=0.5) 
            axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
            legend("bottomleft", legend = round(llik0,3), bty="n", cex=0.7)
            plot(record[1,],ylim=c(0,ylim),type="l", las=2, xaxt="n",main="Newton-Raphson algorithm results", col=my.palette[1],lwd=3, ylab="Value of the variance component", xlab="Iterations processed", cex.axis=0.5) 
            axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
            for(t in 1:(dim(record)[1])){
              lines(record[t,],col=my.palette[t],lwd=3)
            } 
            
            ww <- dim(record)[1]
            lege <- list()
            #lege2 <- list()
            for(k in 1:length(V)){
              if(k == length(V)){
                lege[[k]] <- paste("Var(e):",round(record[k,wi+1],4), sep="")
                #lege2[[k]] <- paste("Var(e):")
              }else{
                lege[[k]] <- paste("Var(",varosss[k],"):",round(record[k,wi+1],4), sep="")
                #lege2[[k]] <- paste("Var(u",k,"):",sep="")
              }
            }
            legend("topleft",bty="n", cex=0.7, col=my.palette, lty=1, lwd=3, legend=unlist(lege))
          }
          ############# END OF LIKELIHOOD PLOT
          
          ## check the change in llik is small
          if(cycle > 1 & abs(delta.llik) < tol*10){
            if(!silent){
              setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
            }
            break
          } 
          if(max(abs(x)) < tol){
            if(!silent){
              setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
            }
            break
          }
        } # END OF ITERATIONS
        
        sigmaxxx <- vector(mode="numeric", length=length(sigma)+length(badd))
        sigmaxxx[badd] <- 0
        sigmaxxx[goodd] <- sigma
        sigma <- sigmaxxx
      }### END OF RECALCULATING VAR.COMP WHEN NEGATIVE
      ###*****************************************************************************************
      
    }else{ ## &&&&&&&&&&&&  "IF FORCED" &&&&&&&&&&&&&&
      
      sigma <- forced
      setTxtProgressBar(pb, (tot/tot))
      
      cat("\nVariance components forced\n")
      var.com2 <- as.matrix(forced, ncol=1)
      logL <- 0
      
      if(is.null(names(ZETA))){
        varosss <- c(paste("u.",df.ord, sep=""))
      }else{
        varosss <- c(names(ZETA))
      }
      #print(var.com2)
    }## ## &&&&&&&&&&&&  END OF FORCING &&&&&&&&&&&&&&
    
    
  }### end of FIXED-ONLY(ONLY INTERCEPTS) vs MIXED model
  
  ################################
  ################################
  # RANDOM VARIABLES
  ### Coeficcient matrix
  ################################
  zvar <- 1:length(ZETA)
  AIC = as.vector((-2 * llik ) + ( 2 * dim(xm)[2]))
  BIC = as.vector((-2 * llik ) + ( log(length(y)) * dim(xm)[2]))
  
  #zvar <- which(unlist(lapply(ZETA, function(x){names(x)[1]})) == "Z")
  
  #####################
  var.com2 <- as.matrix(sigma)
  varo <- as.list(sigma)  # variance components for random, no error
  
  Gspo <- lapply(as.list(c(1:length(ZETA))),function(x,K,v){
    oo=K[[x]]*as.numeric((v[[x]])); return(oo)
  }, K=Gs, v=varo) ##
  Gsp <- as(do.call("adiag1", Gspo),Class="sparseMatrix") # in diagonal
  Rsp <- as(R*as.numeric(var.com2[length(var.com2)]),Class="sparseMatrix")
  
  varo=NULL
  Gspo=NULL
  #########################################################
  # Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467)
  # R-  --  [R-Z[Z'R-Z+G-]-Z'R-]  #-- means minus
  
  vm <- Zsp%*%crossprod(Gsp,tZsp) + Rsp # ZGZ+R
  #vm <- W
  Vinv <- solve(vm,sparse=TRUE, tol = 1e-19)
  #################
  #Vinv2 <- Vinv
  ## Fixed effects
  # B=(X'X V- XVy)-
  xvx <- crossprod(xm, Vinv %*% xm)
  xvxi <- solve(xvx) # variance of fixed effects
  beta <- xvxi %*% crossprod(xm, Vinv %*% y)
  #################
  Var.u <- vector(mode="list", length = length(zvar))
  PEV.u <- Var.u
  u <- Var.u
  
  pm=Vinv-Vinv%*%xm%*%(xvxi%*%txm%*%(Vinv))
  ZETA3 <- lapply(ZETA, function(x){y=list(Z=as(x[[1]],Class="sparseMatrix"),K=as(x[[2]],Class="sparseMatrix"))})
  for(h in zvar){
    Var.u[[h]] <- (as.numeric(var.com2[h,1])^2) *  ( crossprod(ZETA3[[h]][[1]]%*%ZETA3[[h]][[2]], pm)  %*%  (ZETA3[[h]][[1]]%*%ZETA3[[h]][[2]])   ) # sigma^4 ZKP ZK
    PEV.u[[h]] <- as.numeric(var.com2[h,1]) * ZETA3[[h]][[2]] - Var.u[[h]]  # standard errors (SE) for each individual
  } #(sigma2.u^2) *  ( crossprod(Z%*%K, P)  %*%  (Z%*%K)   )
  #################
  ## Random effects
  #u <- list() # we put it up
  ee <-  (y - (xm %*% beta))
  
  if(length(ZETA)==1 & MTG2==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
    ### MTG2 not active currently for NR method
    #for (k in zvar) { # u = KZ'V- (y- XB)
    #  u[[k]] <- solve(t(Us[[k]])) %*% ( ( (ZETA3[[k]][[2]]*as.numeric(var.com2[k,1])) %*% t(ZETA3[[k]][[1]]) %*% Vinv %*% ee ))
    #}
  }else{
    for (k in zvar) { # u = KZ'V- (y- XB)
      u[[k]] <- ( (ZETA3[[k]][[2]]*as.numeric(var.com2[k,1])) %*% t(ZETA3[[k]][[1]]) %*% Vinv %*% ee )
    }
  }
  u <- u[zvar]
  ###############
  #residuals2 <- (y - (xm %*% beta))
  
  
  fitted.u <- 0
  for(h in 1:length(zeta.or)){
    #fitted.y <- fitted.y + (zeta.or[[h]][[1]] %*% u[[h]])
    fitted.u <- fitted.u + (zeta.or[[h]][[1]] %*% u[[h]])
  }
  fitted.y <- (x.or %*% beta) + fitted.u
  fitted.y.good <- fitted.y[good]
  residuals3 <- y - fitted.y[good] # conditional residuals
  ###################
  rownames(beta) <- colnames(xm)
  for(i in 1:length(ZETA)){
    rownames(u[[i]]) <- colnames(ZETA[[i]][[1]])
  }
  ###################
  ####### FISHER's INFORMATION 
  if(Fishers){
    fishers <- A
    fishers.inv <- solve(fishers)
  }else{fishers.inv=NULL}
  ####################
  ###################
  ###################
  if(!is.null(names(ZETA))){
    names(u) <- names(ZETA)
    names(Var.u) <- names(ZETA)
    names(PEV.u) <- names(ZETA)
  }
  logL <- as.vector(llik)
  ###################
  out1 <- as.matrix(var.com2, ncol=1); colnames(out1) <- "Variance Components" # variance components
  rownames(out1) <- c(paste("Var(",varosss,")",sep=""), "Var(Error)")
  res <- list(var.comp=out1, V.inv = Vinv, u.hat=u, Var.u.hat=Var.u, 
              PEV.u.hat=PEV.u, beta.hat=beta, Var.beta.hat=xvxi, 
              LL=logL, AIC=AIC, BIC=BIC, X=xm, fitted.y=fitted.y, 
              fitted.u=fitted.u, residuals=ee, cond.residuals=residuals3,
              fitted.y.good=fitted.y.good, Z=Zsp, K=Ksp, fish.inv=fishers.inv)
  
  layout(matrix(1,1,1))
  #print(logL2.stored)
  #print(big.peaks.col(logL2.stored, -100000))
  return(res)
}


