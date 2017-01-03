###############################################################
###############################################################
## START THE FUNCTION
###############################################################
###############################################################
## define response variable

NRR <- function(y,X=NULL,Z=NULL,R=NULL,tolpar=1e-6,tolparinv=1e-6,maxcyc=10,
                draw=TRUE, constraint=TRUE){
  ### ====== USEFUL FUNCTIONS ====== ###
  '%!in%' <- function(x,y)!('%in%'(x,y))
  NRstep <- function(y,X=NULL,Z=NULL,R=NULL,tolpar=1e-6,tolparinv=1e-6,maxcyc=10,
                     draw=TRUE){
    
    or.var <- var(y,na.rm=TRUE)
    y <- scale(y)
    ## impute y 
    NANA <- which(is.na(y))
    if(length(NANA)>0){y[NANA] <- mean(y,na.rm=TRUE)}
    n <- length(y)#dim(Z1)[1]
    In <- diag(n)
    ## get X
    if(is.null(X)){
      X <- matrix(1,nrow=n)
    }
    ## add the derivatives of the random effects
    if(!is.null(Z)){
      ZKZ <- lapply(Z, function(x){x$Z%*%x$K%*%t(x$Z)})
      ## start the R part
      nz <- length(Z)
    }else{
      ZKZ <- list()
      ## start the R part
      nz <- 0
    }
    ## create derivatives for R matrices
    if(!is.null(R)){
      ttypes <- lapply(R, function(x){x$type[which(x$type!="ID")]})
      nnnn <- lapply(ttypes, function(x){length(which(x!="ID"))}) # no. of R matrices
      nnnn2 <- nnnn
      counter <- 1+nz
      ## keep track of the residual structures
      for(o in 1:length(nnnn)){
        arma.mark <- length(which(ttypes[[o]]=="ARMA"))
        if(length(arma.mark)>0){ # if arma there's an extra parameter lambda to estimate
          nnnn2[[o]] <- seq(counter,counter+nnnn[[o]]+arma.mark-1,1); 
          counter <- counter+nnnn[[o]]+arma.mark
        }else{ # if AR1 or CS
          nnnn2[[o]] <- seq(counter,counter+nnnn[[o]]-1,1); 
          counter <- counter+nnnn[[o]]
        }
      }
      
      ## nnnn tells me how many r structures need derivative
      ## nnnn2 tells me where they whould be stored in ZKZ
      for(i in 1:length(ttypes)){ # for each R structure
        (prot <- ttypes[[i]] )
        (pron <- nnnn2[[i]])
        prono <- pron
        touse <- which(R[[i]]$type != "ID")
        for(j in 1:length(prot)){ # get the derivative for such R structure
          jj <- touse[j]
          if(prot[j]=="AR1"){
            provZKZ <- R[[i]][-length(R[[i]])] # get the R complete except for the "type"
            (nrows <- dim(R[[i]][[jj]])[1])
            (ss1 <- R[[i]][[jj]][1,2])
            # derivative and replace in prov ZKZ which has Rs in order
            provZKZ[[jj]] <- (abs(row(diag(nrows))-col(diag(nrows))) * ss1^(abs(row(diag(nrows)) - col(diag(nrows))) - 1))
            ZKZ[[pron[j]]] <- do.call("kronecker",provZKZ)
          }else if(prot[j]=="CS"){
            provZKZ <- R[[i]][-length(R[[i]])] # get the R complete except for the "type"
            (nrows <- dim(R[[i]][[jj]])[1])
            (ss1 <- R[[i]][[jj]][1,2])
            # derivative and replace in prov ZKZ which has Rs in order
            provZKZ[[jj]] <- matrix(1,nrows,nrows)#(abs(row(diag(nrows))-col(diag(nrows))) * ss1^(abs(row(diag(nrows)) - col(diag(nrows))) - 1))
            ZKZ[[pron[j]]] <- do.call("kronecker",provZKZ)
          }else if(prot[j]=="ARMA"){
            for(rrr in 1:2){ # because there is 2 parameters
              provZKZ <- R[[i]][-length(R[[i]])] # get the R complete except for the "type"
              (nrows <- dim(R[[i]][[j]])[1])
              (ss1 <- R[[i]][[j]][1,2])# ar
              (ss2 <- (R[[i]][[j]][1,3])/ss1) # lam
              # derivative and replace in prov ZKZ which has Rs in order
              amo <- prono[j:(j+1)]
              # derivative for ar
              provZKZ[[j]] <- ARMA.mat.dar(ss1,ss2,nrows)#(abs(row(diag(nrows))-col(diag(nrows))) * ss1^(abs(row(diag(nrows)) - col(diag(nrows))) - 1))
              ZKZ[[amo[rrr]]] <- do.call("kronecker",provZKZ)
            }
            pron <- pron[-j] # delete one so it can atch still
          }#else{
          #  stop(paste("Type specified for the R matrix",j,"in the element",i,"is not available. \nOnly AR1, CS, and ARMA types are available currently."))
          #}
          
        } # end of getting derivative
        
      } # end of derivatives for each R list element
      
      # now ZKZ is complete obtain the R = Ri + ... + Rj
      # where Ri is Rii x Rij , where "x" is kronecker
      Rs <- 0
      for(u in 1:length(R)){
        #hg <- nnnn2[[u]] # matrices for that R located in ZKZ hg
        Ri <- do.call("kronecker",R[[u]][-length(R[[u]])]) # but removing type
        Rs <- Rs + Ri
      }
      ## insert R in the last element of ZKZ
      ZKZ[[length(ZKZ)+1]] <- Rs
    }else{
      ZKZ[[length(ZKZ)+1]] <- In
    }
    
    k <- length(ZKZ) # no. of variance components
    sigma <- rep(var(y),length(ZKZ)) # initial values of variance components
    ## get the initial values for the R structures
    if(!is.null(R)){
      for(kq in 1:length(R)){
        hn <- nnnn2[[kq]]
        ht <- R[[kq]]$type
        ht <- ht[which(ht!="ID")]
        touse <- which(R[[kq]]$type != "ID")
        for(w in 1:length(hn)){
          ji <- hn[w]
          ww <- touse[w]
          if(ht[w]=="AR1"){
            sigma[ji] <- R[[kq]][[ww]][1,2] # this is the initial var.comp
          }else if(ht[w]=="CS"){
            sigma[ji] <- R[[kq]][[ww]][1,2] # this is the initial var.comp
          }else if(ht[w]=="ARMA"){
            (sigma[ji] <- R[[kq]][[ww]][1,2])# ar
            (sigma[ji+1] <- (R[[kq]][[ww]][1,3])/sigma[ji]) # lam
            ht <- ht[-ww]
          }
        }
      }
    }
    
    coef <- sigma
    pos <- rep(FALSE,k)
    
    qr <- qr(X)
    rankQ <- n-qr$rank
    X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
    rankQX <- n - qr$rank
    logL2.stored <- numeric()
    
    taper <- rep(0.9, maxcyc) # weighting parameter
    taper[1:2] <- c(0.5, 0.7)
    
    
    ## ==============
    ####### algorithm
    ## ==============
    #Rmarka <- FALSE
    #NOW <- FALSE
    sigma.old <- as.matrix(sigma)
    #maxcyc <- 40
    
    for(cycle in 1:maxcyc){
      #V = ZKZ + ZKZ + ... + ZRZ
      
      Sigma <- 0
      if(!is.null(Z)){
        hhh <- c(1:nz,length(sigma))
      }else{
        hhh <- c(length(sigma))
      }
      
      for(i in hhh) { # 1:length(sigma)
        Sigma <- Sigma + ZKZ[[i]]*sigma[i]
      }
      
      ### invert V
      V <- try(solve(as(Sigma, Class="sparseMatrix"),In,sparse=TRUE), silent = TRUE)
      if(class(V) == "try-error"){
        V <- try(solve(as(Sigma + tolparinv * diag(dim(Sigma)[2]), Class="sparseMatrix"),In,sparse=TRUE), silent = TRUE)
      }#V <- solve(Sigma,In)
      #print(V[1:5,1:5])
      # P = V-  -  V-X [X' V- X]-1 XV- =  WQK
      VX <- V %*% X # V- X
      P <- V - VX %*% solve(t(X)%*%VX, t(VX)) #WQK
      # y'Py
      rss <- as.numeric(t(y) %*% P %*% y)
      
      #rssR <- as.numeric(t(y) %*% PR %*% y)
      # pos is rep(FALSE,k)
      sigma <- sigma * rss/rankQX
      coef[!pos] <- sigma[!pos] # FALSE are copied
      coef[pos] <- log(sigma[pos]) # TRUE are copied
      P <- P * rankQX/rss # P * [r(X) / yPy] #WQX <- WQX * rankQK/rss
      rss <- rankQX 
      eig <- sort(eigen(P,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankQX]
      
      if(any(eig < 0)){
        P <- P + (tolpar - min(eig))*diag(dim(P)[1])
        eig <- eig + tolpar - min(eig)
      }
      
      ldet <- sum(log(eig)) 
      llik <- ldet/2 - rss/2 #.5 [log(det(eigen(P))) - yPy]
      if(cycle == 1) llik0 <- llik #keep first log likelihood
      delta.llik <- llik - llik0 # increase of likelihood with respect to initial LL
      llik0 <- llik # update likelihood to current iteration
      
      logL2.stored[cycle] <- llik
      
      x <- NULL # a clean x
      var.components <- rep(1,k) # variance components
      ind <- which(pos) # which are TRUE
      if(length(ind)) var.components[ind] <- sigma[ind] #update
      
      TT <- list() #list of PVi  Vi=dV/ds
      ## P Vi
      TT[[k]] <- P %*% ZKZ[[k]]# last element (for error) only P
      if(k>1) { #which is the case, must be the number corresponding to error
        # i.e. 3:1   P ZKZ = P Vi
        for(ii in (k-1):1) TT[[ii]] <- P %*% ZKZ[[ii]] ## # P ZKZ' #... V.inv - V.inv X [X' V.inv X]-1 * ZKZi
      }
      
      ## obtain first derivatives
      # Vi = dV/ds  ..... y' P Vi P y - tr(P Vi)  same than -tr(PVi) - y'PViPy
      x <- sapply(TT,function(x) as.numeric(t(y) %*% x %*% P %*% y - sum(diag(x))))
      ## theta(k) * dL/ds  ..... are scalar values
      x <- x * var.components
      ## second derivatives .... [theta(i) * 1st.deriv(i)] * [theta(j) * 1st.deriv(j)]  * sigma(i) * sigma(j)
      A <- matrix(rep(0, k^2), k, k)
      entries <- expand.grid(1:k,1:k) # indices to be filled
      ## Fisher's Information tr(PA*PA*) .... A*=Vi=dV/ds .... [Vi Vj'] si sj
      #print(var.components)
      ff <- function(x) sum(TT[[x[1]]] * t(TT[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
      aa <- apply(entries,1,ff) # matrix of combinations of var.comp 
      A[as.matrix(entries)] <- aa
      ######### one including Rs other excluding them
      A.svd <- ginv(A) # Inverse of Fishers
      x <- A.svd %*% x #update 
      coef <- coef + taper[cycle] * x # sigma + f[s*F-*dL/ds] ..... = coef + taper[x]
      sigma[!pos] <- coef[!pos] # FALSES are replaced
      sigma[pos] <- exp(coef[pos]) # TRUES are replaced
      
      #print(sigma)
      # for the first 7 cycles don't update Ri's
      #if(cycle <= 0){
      #  sigma[3:4] <- sigma.old[3:4]
      #  #sigma[4] <- sigma.old[4]
      #}
      
      if(cycle > 1 & abs(delta.llik) < tolpar) {
        break
      }
      
      sigma.old <- cbind(sigma.old,as.matrix(sigma))
      
      ############ =================================== #############
      ############ =================================== #############
      ############ =================================== #############
      ## update R and ZKZ(derivatives for Rs)
      if(!is.null(R)){
        for(qq in 1:length(nnnn2)){ # for each R
          mii <- nnnn2[[qq]]
          sisi <- sigma[mii]
          # now update R
          touse <- which(R[[qq]]$type != "ID")
          for(p in 1:length(sisi)){ # for each Rii
            pp <- touse[p]
            momo <- R[[qq]]$type[pp] # is AR1, CS ARMA??
            ## possibilities
            if(momo == "AR1"){
              ## update Rii
              R[[qq]][[pp]] <- AR1.mat(sisi[p],dim(R[[qq]][[pp]])[1])
            }else if(momo=="CS"){
              R[[qq]][[pp]] <- CS.mat(sisi[p],dim(R[[qq]][[pp]])[1])
            }else if(momo=="ARMA"){
              
            }else if(momo=="ID"){
              R[[qq]][[pp]] <- ID.mat(1,dim(R[[qq]][[pp]])[1])
            }
            ##################
            
          }
          ### end of update Ri
        }
        
        ## end of update Rs
        ############ =================================== #############
        ############ =================================== #############
        ############ =================================== #############
        #                  NOW UPDATE DERIVATIVES
        for(qq in 1:length(nnnn2)){
          mii <- nnnn2[[qq]]
          sisi <- sigma[mii]
          # now update R
          touse <- which(R[[qq]]$type != "ID")
          for(p in 1:length(sisi)){
            pp <- touse[p]
            momo <- R[[qq]]$type[pp] # is AR1, CS ARMA??
            ## possibilities
            if(momo == "AR1"){
              ## update Rii
              ## update its derivative
              provZKZ <- R[[qq]][-length(R[[qq]])] # get the R complete except for the "type"
              (nrows <- dim(R[[qq]][[pp]])[1])
              (ss1 <- sisi[p])
              ko <- mii[p]#nnnn2[[q]][p]
              # derivative and replace in prov ZKZ which has Rs in order
              provZKZ[[pp]] <- sigma[length(sigma)]*(abs(row(diag(nrows))-col(diag(nrows))) * ss1^(abs(row(diag(nrows)) - col(diag(nrows))) - 1))
              ZKZ[[ko]] <- do.call("kronecker",provZKZ)
            }else if(momo=="CS"){
              ## update Rii
              ## update its derivative
              provZKZ <- R[[qq]][-length(R[[qq]])] # get the R complete except for the "type"
              (nrows <- dim(R[[qq]][[pp]])[1])
              (ss1 <- sisi[p])
              ko <- mii[p]#nnnn2[[q]][p]
              # derivative and replace in prov ZKZ which has Rs in order
              provZKZ[[pp]] <- sigma[length(sigma)]*matrix(1,nrows,nrows)#(abs(row(diag(nrows))-col(diag(nrows))) * ss1^(abs(row(diag(nrows)) - col(diag(nrows))) - 1))
              ZKZ[[ko]] <- do.call("kronecker",provZKZ)
            }else if(momo=="ARMA"){
              
            }
            ##################
            
          }
          ### end of update Ri
        }
        ############ =================================== #############
        ############ =================================== #############
        ############ =================================== #############
        ## update REAL R after updating Ri's and derivatives
        Rs <- 0
        for(uu in 1:length(R)){
          #hg <- nnnn2[[u]] # matrices for that R located in ZKZ hg
          Ri <- do.call("kronecker",R[[uu]][-length(R[[uu]])]) # but removing type
          Rs <- Rs + Ri
        }
        ## insert R in the last element of ZKZ
        ZKZ[[length(ZKZ)]] <- Rs
      }# END OF if(is.null(R))
      ##########################
      if(draw){
        plot(logL2.stored,type="o",xlab="Iterations",col="cadetblue")
      }
    }# end of iterations for variance components
    
    ##
    if(!is.null(R)){
      rstr <- unlist(nnnn2)
      zstr <- setdiff(1:length(sigma), rstr)
      sigma2 <- sigma
      sigma2[zstr] <- sigma2[zstr]*or.var
    }else{
      rstr <- NULL
      zstr <- 1:length(sigma)
      sigma2 <- sigma
      sigma2[zstr] <- sigma2[zstr]*or.var
    }
    
    return(list(sigma2,rstr,zstr,llik,A,A.svd,ZKZ,sigma,pos))
  }
  
  #print(str(Z))
  #print(str(R))
  ### ====== calculate var.comp ====== ###
  sig <- NRstep(y=y,X=X,Z=Z,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=maxcyc,
                draw=draw)
  #print(sig[[1]])
  ### ====== constraint for zero var.comp ====== ###
  if(constraint){
    varo <- sig[[3]] # this are var.comp G which should be positive
    bad <- which(sig[[1]][varo] < 0)
    if(length(bad)>0){ #there's negative variance components which should be zero
      realbad <- varo[bad]
      # error variance was negative and is the only one negative
      # don't do anything
      if( (c(length(sig[[1]])) %in% realbad) & (length(realbad)==1)){ 
        sig3 <- sig
        sigi <- sig[[1]]
      }else if((length(sig[[1]]) %in% realbad) & (length(realbad)>1)){
        # error was negative but along with others, only recaluclate the others
        realbad2 <- realbad[-length(realbad)] # without R or sigma.e
        if(realbad2 == length(Z)){ # all random effects are zero
          ZZ <- NULL
        }else{ #only some random effects are zero
          ZZ <- Z[-realbad2]
        }
        good <- setdiff(1:length(sig[[1]]),realbad2)
        sig2 <- NRstep(y=y,X=X,Z=ZZ,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=maxcyc,
                       draw=draw)
        sigo <- sig[[1]]
        sigo[good] <- sig2[[1]]
        sigo[realbad2] <- 0
        sig3[[1]] <- sigo
      }else if((length(sig[[1]]) %!in% realbad) & (length(realbad)>0)){
        # error was not negative but, only re-calculate the others
        sig3 <- sig
        realbad2 <- realbad # only Gs
        if(length(realbad2) == length(Z)){ # all random effects are zero
          ZZ <- NULL
        }else{ #only some random effects are zero
          ZZ <- Z[-realbad2]
        }
        good <- setdiff(1:length(sig[[1]]),realbad2)
        sig2 <- NRstep(y=y,X=X,Z=ZZ,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=maxcyc,
                       draw=draw)
        sigo <- sig[[1]]
        sigo[good] <- sig2[[1]]
        sigo[realbad2] <- 0
        sig3[[1]] <- sigo
      }
    }else{ # if constraint but no negative var.comp
      sig3<-sig
    }
  }else{
    sig3<-sig
  } # if no constraint
  ## sig3 has 1) var.comp, 2) Ri structures, 3) G struct and sigma.e, 4) llik, 5) F, 6) F.inv
  ## 7) ZKZ, 8) sigma scaled, 9) pos
  
  ### ====== get other parameters ====== ###
  n <- length(y)
  NANA <- which(is.na(y));if(length(NANA)>0){y[NANA] <- mean(y,na.rm=TRUE)}
  if(is.null(X)){
    xm <- matrix(1,nrow=n)
  }else{
    xm <- as.matrix(X) 
  };txm <- t(xm)
  
  zvar <- sig3[[3]]; zvar <- zvar[-length(zvar)]
  var.com2 <- as.matrix(sig3[[1]])
  llik <- sig3[[4]]
  hhh <- sig3[[3]]
  ZKZ <- sig3[[7]]; #sig3[[7]]<- NULL
  
  ################################
  ################################
  # RANDOM VARIABLES
  ### Coeficcient matrix
  ################################
  AIC = as.vector((-2 * llik ) + ( 2 * dim(xm)[2]))
  BIC = as.vector((-2 * llik ) + ( log(length(y)) * dim(xm)[2]))
  #####################
  
  Sigma <- 0
  for(i in hhh) { # 1:length(sigma)
    Sigma <- Sigma + ZKZ[[i]]*var.com2[i,1]
  }
  Vinv <- try(solve(as(Sigma, Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
  if(class(Vinv) == "try-error"){
    Vinv <- try(solve(as(Sigma + tolparinv * diag(dim(Sigma)[2]), Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
  }
  
  #################
  ## Fixed effects
  # B=(X'X V- XVy)-
  xvx <- crossprod(xm, Vinv %*% xm)
  xvxi <- solve(xvx) # variance of fixed effects
  beta <- xvxi %*% crossprod(xm, Vinv %*% y)
  #################
  Var.u <- vector(mode="list", length = length(zvar))
  PEV.u <- Var.u
  u <- Var.u
  ## fill the lists
  pm=Vinv-Vinv%*%xm%*%(xvxi%*%txm%*%(Vinv))
  if(!is.null(Z)){
    Z <- lapply(Z, function(x){y=list(Z=as(x[[1]],Class="sparseMatrix"),K=as(x[[2]],Class="sparseMatrix"))})
    
    for(h in zvar){
      Var.u[[h]] <- (as.numeric(var.com2[h,1])^2) *  ( crossprod(Z[[h]][[1]]%*%Z[[h]][[2]], pm)  %*%  (Z[[h]][[1]]%*%Z[[h]][[2]])   ) # sigma^4 ZKP ZK
      PEV.u[[h]] <- as.numeric(var.com2[h,1]) * Z[[h]][[2]] - Var.u[[h]]  # standard errors (SE) for each individual
    } #(sigma2.u^2) *  ( crossprod(Z%*%K, P)  %*%  (Z%*%K)   )
  }
  
  #################
  ## Random effects
  #u <- list() # we put it up
  ee <-  (y - (xm %*% beta))
  
  if(!is.null(Z)){
    for (aaa in zvar) { # u = KZ'V- (y- XB)
      u[[aaa]] <- ( (Z[[aaa]][[2]]*as.numeric(var.com2[aaa,1])) %*% t(Z[[aaa]][[1]]) %*% Vinv %*% ee )
    }
    u <- u[zvar]
  }
  ###############
  #residuals2 <- (y - (xm %*% beta))
  
  
  fitted.u <- 0
  if(!is.null(Z)){
    for(h in 1:length(Z)){
      #fitted.y <- fitted.y + (zeta.or[[h]][[1]] %*% u[[h]])
      fitted.u <- fitted.u + (Z[[h]][[1]] %*% u[[h]])
      ## in other things
      rownames(u[[h]]) <- colnames(Z[[h]][[1]])
    }
  }
  fitted.y <- (xm %*% beta) + fitted.u
  residuals3 <- y - fitted.y # conditional residuals
  ###################
  rownames(beta) <- colnames(xm)
  ###################
  ####### FISHER's INFORMATION 
  ## sig3 has 1) var.comp, 2) Ri structures, 3) G struct and R, 4) llik, 5) F, 6) F.inv
  ## 7) ZKZ
  A.svd <- sig3[[6]]
  A <- sig3[[5]]
  sigma <- sig3[[8]]
  pos <- sig3[[9]]
  
  sigma.cov <- (A.svd * 2) 
  FI <- A/2
  ## convert FI using pos
  FI.c <- matrix(0,dim(FI)[1],dim(FI)[2])
  FI.c <- FI / tcrossprod((sigma-1)*pos+1)
  #print(A.svd)
  #names(sigma) <- Vcoef.names
  sigma.cov <- try(ginv(FI.c),silent=TRUE)
  error1 <- (class(sigma.cov)=="try-error")
  #if(error1) {
  #  cat("Warning: solution lies on the boundary; check sigma & pos\nNo standard errors for variance components returned\n")
  #  sigma.cov <- matrix(NA,k,k)
  #}
  
  ###################
  if(!is.null(names(Z))){
    names(u) <- names(Z)
    names(Var.u) <- names(Z)
    names(PEV.u) <- names(Z)
  }
  logL <- as.vector(llik)
  ###################
  
  #if(is.null(forced)){
  # rownames(sigma.cov) <- colnames(sigma.cov) <- Vcoef.names
  #}
  
  ## Additional warning if any pos entries are TRUE and the corresponding term is close to zero
  ## Only give this warning it I haven't given the previous one.
  
  
  #out1 <- as.matrix(var.com2, ncol=1); colnames(out1) <- "Variance Components" # variance components
  #rownames(out1) <- c(paste("Var(",varosss,")",sep=""), "Var(Error)")
  
  #if(is.null(forced)){
  #  rownames(sigma.cov) <- rownames(out1)
  #  colnames(sigma.cov) <- rownames(out1)
  #}
  if(is.null(rownames(var.com2))){
    rownames(var.com2) <- paste("u",1:length(var.com2),sep=".")
    ## add names to RE
    if(!is.null(Z)){ # there's RE
      if(is.null(names(Z))){
        rownames(var.com2)[zvar] <- paste("u",1:length(Z),sep=".")
      }else{
        rownames(var.com2)[zvar] <- names(Z)
      }
    }
    #############
    ## add names to r structures
    if(!is.null(R)){ # thre's r structures
      rvaro <- sig3[[2]]
      rownames(var.com2)[rvaro] <- paste("R.cor",1:length(rvaro),sep=".")
    }
    ##############
    ## add sig,a.e names
    rownames(var.com2)[length(var.com2)] <- "Residual"
  }
  
  res <- list(var.comp=var.com2, V.inv = Vinv, u.hat=u, Var.u.hat=Var.u, 
              PEV.u.hat=PEV.u, beta.hat=beta, Var.beta.hat=xvxi, 
              LL=logL, AIC=AIC, BIC=BIC, X=xm, fitted.y=fitted.y, 
              fitted.u=fitted.u, residuals=ee, cond.residuals=residuals3,
              fish.inv=sigma.cov, sigma.scaled=sigma, forced=FALSE)
  
  return(res)
}

