AI <- function(y, X=NULL, ZETA=NULL, R=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=50, constraint=TRUE, init=NULL, sherman=FALSE, che=TRUE, EIGEND=FALSE, Fishers=FALSE, gss=TRUE, forced=NULL, tolpar = 1e-04, tolparinv = 1e-06){
  
  convergence <- FALSE
  if(EIGEND){
    DISO <- dim(ZETA[[1]]$Z)
    if(DISO[1] != DISO[2]){
      stop("EIGEN DECOMPOSITION EIGEND ONLY WORKS FOR SQUARE PROBLEMS 
           'Z' MATRIX IS NOT SQUARE",call.=FALSE)
      
    }
  }
  
  y.or <- y
  x.or <- X
  pushess <- (rep(0,length(ZETA)+1))
  ### make full function
  '%!in%' <- function(x,y)!('%in%'(x,y))
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  ### zig.zag function
  zig.zag <- function(zz){
    prove <- diff(zz)
    res <- vector()
    for(i in 2:length(prove)){
      if((prove[i] < 0 & prove[i-1] > 0) | (prove[i] > 0 & prove[i-1] < 0)){
        res[i] <- TRUE
      }else{
        res[i] <- FALSE
      }
    }
    res <- res[-1]
    probab <- length(which(res))/length(res)
    return(probab)
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
    yv <- scale(y)
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
    if(is.null(R)){R <- list(units=diag(length(y)))} # dim(x[[1]])[2]
    
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
          }
        }else{y <- x}; 
        return(y)
      })
    }
    #######################################################
    ## order random effects according to degrees of freedom
    tokeep <- names(ZETA)
    #     df <- unlist(lapply(ZETA, function(x){dim(x[[1]])[2]}))
    #     df2 <- sort(df, decreasing = FALSE)
    #     df.ord <- numeric() # # 
    #     for(u in 1:length(df)){
    #       df.ord[u] <- which(df2 %in% df[u])[1]
    #       df2[df.ord[u]] <- NA
    #     }
    #     ZETA <- ZETA[df.ord]
    #     names(ZETA) <- tokeep[df.ord]
    if(!is.null(init)){init <- init}#[c(df.ord, (length(init)))]}
    if(!is.null(forced)){forced <- forced}#[c(df.ord, (length(forced)))]}
    #####################################################
    ## to use later for fitted values
    x.or <- as.matrix(xm)
    zeta.or <- ZETA
    zeta.or  <- lapply(zeta.or , function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    R.or <- R
    ##
    if((length(ZETA)==1) & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2]) & EIGEND){
      misso <- which(is.na(y))
      if(length(misso) >0){
        y[misso] <- median(y, na.rm=TRUE)
      }
    }
    ZETA2 <- ZETA; y2 <- y ; good <- which(!is.na(y)) # make a copy and identify no missing data
    R2 <- R
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
      
      R <- lapply(R2, function(x,good){x <- x[good,good]; return(x)}, good=good)
    }else{
      ZETA <- lapply(ZETA2, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
      R <- lapply(R2, function(x,good){x <- x[good,good]; return(x)}, good=good)
    }
    ################
    y <- y[good]
    ZETA <- lapply(ZETA, function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    xm <- as.matrix(xm[good,])
    txm <- t(xm)
    #R <- R[good,good]
    qr <- qr(xm)
    rankX <- length(good)-qr$rank
    ###########################
    #### BECOME SPARSE
    #R=as(R,Class="sparseMatrix")
    if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
      EIGENS <- lapply(ZETA, function(x){eigen(x[[2]])}) # eigen decomposition of K
      Us <- lapply(EIGENS, function(x){x$vectors}) # extract eigen vectors U
      Usp <- as(do.call("adiag1", Us),Class="sparseMatrix") # U'G as diagonal
      Ds <- lapply(EIGENS, function(x){diag(x$values)}) # extract eigen values D
      Dsp <- as(do.call("adiag1", Ds),Class="sparseMatrix") # U'G as diagonal
      ZETA <- lapply(as.list(1:length(ZETA)),function(x,zz,kk){list(Z=zz[[x]][[1]], K=kk[[x]])}, zz=ZETA, kk=Ds)
    }
    Zs <- lapply(ZETA, function(x){x[[1]]})
    Gs <- lapply(ZETA, function(x){x[[2]]})
    Zsp <- as(do.call("cbind", Zs),Class="sparseMatrix") # column bind Z=[Z1 Z2 Z3]
    tZsp <- t(Zsp)
    Ksp <- as(do.call("adiag1", Gs),Class="sparseMatrix") # G as diagonal
    fail=FALSE
    ##########################################################
    # We fill fill a list with derivatived to be used in the algorithm
    # om has ZKZ elements which are the derivatives of V with respect to the variance comp.
    # dV/ds = d(ZKsZ + Ie)/d(s) = ZKZ
    ZETA2 <- lapply(ZETA, function(x){y=list(Z=as(x[[1]],Class="sparseMatrix"),K=as(x[[2]],Class="sparseMatrix"))})
    zvar <- which(unlist(lapply(ZETA, function(x){names(x)[1]})) == "Z")
    om <- list()
    for (k in zvar) {
      om[[k]] <- tcrossprod(ZETA2[[k]][[1]], ZETA2[[k]][[1]] %*% (ZETA2[[k]][[2]]) ) 
    }
    for(k1 in 1:length(R)){ # change
      om[[k1 + length(zvar)]] <- as(R[[k1]],Class="sparseMatrix") 
    }
    #######################
    ## Initial values
    if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
      y <- as.vector(t(Usp) %*% as.matrix(y,ncol=1))
      xm <- t(Usp) %*% xm
      txm <- t(xm)
    }
    var.y <- var(y, na.rm=TRUE)
    yv <- scale(y)
    nvarcom <- length(ZETA) + length(R) #modified 2.8
    base.var <- var(yv, na.rm = TRUE)#/nvarcom
    ### EIGEND requires very small initial var.comp compared to regular models
    if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[1] == dim(ZETA[[1]][[2]])[2])){
      if(is.null(init)){var.com <- c(rep(.0001, nvarcom))}else{var.com <- init/var.y} # at the end error variance
      #print(var.com)
    }else{ # before instead of base.var was .01
      if(is.null(init)){
        var.com <- c(rep(base.var, nvarcom))
        #provide different initial values to residuals
        #var.com[setdiff(1:length(var.com),1:length(ZETA))] <- var.com[setdiff(1:length(var.com),1:length(ZETA))]*.9
      }else{
        var.com <- init/var.y
      } # at the end error variance
    }
    weird=FALSE
    tn = length(yv)
    logL2=-10000000 # initial log-likelihood
    logL2.stored <- round(logL2,0)
    conv=0 # convergence
    wi=0 # counter
    record<- matrix(var.com*var.y, ncol=1)
    taper <- rep(0.9, iters) # weighting parameter for updates
    taper[1:2] <- c(0.5, 0.7) # c(0.5, 0.7)
    #var.com <- rep(var(yy2, na.rm = TRUE)/nvarcom, nvarcom) ### NEW at the end error variance
    ##################
    if(is.null(names(ZETA))){
      varosss <- c(paste("u.",1:length(ZETA), sep=""))
    }else{
      varosss <- c(names(ZETA))
    }
    ### ADDITION FOR sommer 2.8
    if(is.null(names(R))){
      varosss <- c(varosss,paste("Res.",1:length(R),sep=""))
    }else{
      varosss <- c(varosss,names(R))
    }
    
    lege2 <- list()
    for(k in 1:length(var.com)){ # 
      gh1 <- varosss[k]
      if(k == length(var.com)){
        lege2[[k]] <- paste("Var(Residual):")
      }else{
        lege2[[k]] <- paste("Var(",gh1,"):",sep="")
      }
    }
    ##################
    ## initialize the progress bar
    if(!silent){
      count <- 0
      tot <- 15
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    #####################
   
    #########################################################################################################
    ## this portion checks if there is random effects with n < p (less observations than parameters to estimate)
    ## to take the right values in AI when bad values are found
    nmbb <- TRUE
    dimen.zeta <- unlist(lapply(ZETA, function(x){lapply(x,function(y){d1<-dim(y);re<-0;if(d1[2]>d1[1]){re<-1};return(re)})}))
    mbb.check <- which(dimen.zeta > 0)
    if(length(mbb.check)>0){ # if n < p exist
      nmbb <- FALSE
      ### if n < p means very small variance components, so we adjust the initial var.comp
      if(is.null(init)){var.com <- c(rep(.0001, nvarcom))}else{var.com <- init/var.y} # at the end error variance
    }; #print(nmbb)
    #########################################################################################################
    
    if(!is.null(forced)){
      iters<-2 
    }
    ups <- numeric()
    #if(is.null(forced)){ # &&&&&&&&&&&&&&&& IF NOT FORCED &&&&&&&&&&&&&&&&&
    
    while (conv==0) { # ==================== START AI ALGORITHM =========================
      wi=wi+1
      if(!silent){
        count <- count + 1
      }
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # V matrix (page 290 in AJHG 96:283-294) actually G matrix in MME
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      varo <- as.list(var.com)  # variance components as list, no error included
      Gspo <- lapply(as.list(c(1:length(ZETA))),function(x,K,v){
        oo=K[[x]]*as.numeric((v[[x]])); return(oo)
      }, K=Gs, v=varo) ## K*v(u)
      Gsp <- as(do.call("adiag1", Gspo),Class="sparseMatrix") # G as diagonal
      
      Rsp <- R[[1]] * as.numeric(var.com[1 + length(zvar)]) # change
      if(length(R)>1){
        for(tu in 2:length(R)){
          Rsp <- Rsp + as(R[[tu]] *as.numeric(var.com[tu + length(zvar)]),Class="sparseMatrix")
        }
      }
      
      varo <- NULL
      Gspo <- NULL
      #########################################################
      # Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467)
      # R-  --  [R-Z[Z'R-Z+G-]-Z'R-]  #-- means minus
      if(sherman){
        Rinv=solve(Rsp, sparse=TRUE, tol = 1e-19)
        Ginv=solve(Gsp, sparse=TRUE, tol = 1e-19)
        ZRZG= solve( as(tZsp%*%Rinv%*%Zsp + Ginv, Class="sparseMatrix"),sparse=TRUE, tol = 1e-19  )
        vm <- Zsp%*%(Gsp%*%tZsp) + Rsp # V=ZGZ+R
        vmi = Rinv - ( Rinv%*%Zsp%*%ZRZG%*%t(Zsp)%*%Rinv)
        #########################################################
      }else{
        vm <- Zsp%*%crossprod(Gsp,tZsp) + Rsp # V=ZGZ+R, was Gsp %*%t(Zsp)
        vmi <- solve(vm, sparse=TRUE, tol = tolparinv) # inverse of V
      }
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # P matrix (page 290 in AJHG 96:283-294), projection matrix
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      xvx <- crossprod(xm, vmi %*% xm)
      pm <- vmi - vmi %*% xm %*% solve(xvx, crossprod(xm, vmi))
      ytPy <- t(yv)%*%(pm%*%yv)
      
      ## regress tricks!!!!!!!!!
      var.com <- var.com * (as.numeric(ytPy)/as.numeric(rankX))
      pm <- pm * (as.numeric(rankX)/as.numeric(ytPy))

      vmi=NULL
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # log Likelihood (page 290 in AJHG 96:283-294)
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # covariance matrices need to be full rank or negative values will indicate matrices that are not positive semidefinite
      if(REML){ # ======= ********* WHEN REML ************ ==================
        ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
        if(is.infinite(ddv)){ # ======= IF DETERMINANT IS INFINITE(REML) =======
          stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.",call.=FALSE)
        }else{ # ======= IF EVERYTHING GOES WELL =========
          logL=as.numeric(-0.5*((ddv)+determinant(xvx, logarithm = TRUE)$modulus[[1]]+ytPy)) # log likelihood, problem
        }
      }else{ # ======= ********* WHEN ML ************ ==================
        ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
        if(is.infinite(ddv)){ # ======= IF DETERMINANT IS INFINITE (ML) =======
          stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.",call.=FALSE)
        }else{ # ======= IF EVERYTHING GOES WELL(ML) =========
          logL=as.numeric(-0.5*((ddv) + ytPy)) # log likelihood, problem
        }
      }
      #################################
      ### CONTROL OF ZIG.ZAG LIKELIHOOD

      #################################
      
      #############
      ############# WAS abs(logL-logL2)<0.001     
      if (abs(logL-logL2) < tolpar | wi == iters ){ ## CONVERGENCE, not sure if should be absolute value or not
        conv=1
        
        if (abs(logL-logL2) < tolpar){
          convergence <- TRUE
        }
        if(!silent){
          setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
        }
        if(wi == iters){ # LAST RESOURCE CHANGE TO EM ALGORITHM
          boundary <- which(var.com<=0)
        }
      }else{
        if(!silent){
          setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
        }
        logL2=(logL) # replace the initial logLik value for the value reached
        logL2.stored <- c(logL2.stored,logL2)
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        # Fill AI matrix (page 290 in AJHG 96:283-294 or Eq.8 in GSE 38: 25-43)
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        # average of second derivatives and expectations
        # om has ZKZ elements which are the derivatives of V with respect to the variance comp.
        # dV/ds = d(ZKsZ + Ie)/d(s) = ZKZ, this are called Hi terms in Gilmour et al. (1995)
        aim=matrix(0,nvarcom,nvarcom) # average information matrix
        py=pm%*%yv # Py
        
        for(i in 1:nvarcom){ ### (A9) equation in Lee (2015) 
          for(j in 1:i){
            aim[i,j]=as.numeric(0.5*(t(yv)%*%om[[i]]%*%pm%*%om[[j]]%*%pm%*%py)) # .5 * y' Hi %*% P %*% Hj %*% P %*% P %*% y
            # where Hi=ZKZ, here om[[i]]
            if(i != j){
              aim[j,i] <- aim[i,j]
            }
          }
        }
        #print(aim)
        
        if(rankMatrix(aim)[1] == dim(aim)[2]){
          aimi=solve(aim) # solve?? (AI)- # Average Information Matrix Inverse (aimi)
        }else{
          aimi=ginv(aim) # solve?? (AI)- # Average Information Matrix Inverse (aimi)
        }
        
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        # D = matrix of first derivatives (((dL/ds Eq. A10 in Lee (2015) )))
        # (page 290 in AJHG 96:283-294 or Eq.9 in GSE 38: 25-43)
        # update values of Variance Components
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        dldv=matrix(0,nvarcom) # matrix of first derivatives, single column, scores to update variance components
        ##########
        for(k in 1:nvarcom){
          prm=pm%*%om[[k]] # PHi# ( P(e) = Vinv - [Vinv X (X'V-X)- X Vinv]  )  %*% I or K, etc       from I sigma(e) initial matrix for errors
          #tr1=0
          #for (i in 1:tn) { # trace P(e)
          #  tr1=tr1+prm[i,i] # add diagonal values from P(e), trace 1
          #}
          tr1=sum(diag(prm)) # trace
          # fill the first derivative matrix
          dldv[k,1]=-0.5*tr1+0.5*as.numeric(t(yv)%*%prm%*%py) # -.5 * [ tr(PHi)+0.5 * y' %*% PHi %*% Py ]        ## fill 1st derivative matrix 1,1
        }
        #print(aim)
        #print(dldv)
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        # AI update (page 290 in AJHG 96:283-294 or Eq.7 in GSE 38: 25-43)
        
        up=aimi%*%dldv
        #print(up)
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        
        
        vcomp <- var.com
        var.com <- as.matrix(var.com) + (taper[wi]*as.matrix(up))# taper is a regress trick
        
        
        ups <- cbind(ups,var.com)
        #print(var.com)
        # CONSTRAINTS: variance components cannot be zero or exceed total variance
        # if after 5 iterations keeps pushing to negative values make it zero
        ####***********************************
        ####***********************************
        ######### modified for sommer 2.7
        #### the idea is that vc that are negative are set to almost zero
        #### but negative  to make sure that is refitted
        #print(var.com)
        
        if(constraint){
          #print("yes")
          nmm <- length(var.com)
          boundary <- which(var.com<=0)
          var.com[which(var.com<=0),] <- mean(var.com[(nmm-length(R)+1):nmm]) * (1.011929e-07^2)# -(1e-6) #vcomp[which(var.com<=0)]#
          
        }
        ####***********************************
        ####***********************************
        
        record <- cbind(record, var.com * as.numeric(var.y))
        
        
        ###########################################
        ###########
        # PLOT
        if(draw){
          ylim <- max(unlist(record), na.rm=TRUE)
          my.palette <- brewer.pal(7,"Accent")
          layout(matrix(1:2,2,1))
          plot(logL2.stored[-1],type="l", main="logLikelihood", col=my.palette[7],lwd=3, las=2, xaxt="n", ylab="logLikelihood value", xlab="Iterations processed", cex.axis=0.5) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          legend("bottomleft", legend = round(logL2,3), bty="n", cex=0.7)
          plot(record[1,],ylim=c(0,ylim),type="l", las=2, xaxt="n",main="Average Information algorithm results", col=my.palette[1],lwd=3, ylab="Value of the variance component", xlab="Iterations processed", cex.axis=0.5) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          for(t in 1:(dim(record)[1])){
            lines(record[t,],col=my.palette[t],lwd=3)
          } 
          
          ww <- dim(record)[1]
          lege <- list()
          for(k in 1:length(var.com)){
            lege[[k]] <- paste("Var(",varosss[k],"):",round(record[k,wi+1],4), sep="")
          }
          legend("topleft",bty="n", cex=0.7, col=my.palette, lty=1, lwd=3, legend=unlist(lege))
        }
        
        ##$$$$$$$$$$$$$$$$

        #print(abnormal)
        ##$$$$$$$$$$$$$$$$
        
        ###########
        fail=FALSE
        ###########
      }
      
    }# =====================================  END =======================================
    #####################################################################################
    ####### FOR LOOP FOR ITERATIONS FINISHED ############################################
    
    #boundary <- which(var.com==mean(var.com[(nmm-length(R)+1):nmm]) * (1.011929e-07^2))
    if(!is.null(forced)){
      var.com2 <- as.matrix(forced)
    }else{
      var.com2 <- var.com*var(y.or,na.rm=TRUE)
    }
    #########################################################
    # perform boosting if there was a var comp close to zero
    # and AI is able to converge
    
    
    
    ###### ENF OF APPLY CONSTRAINT 
    ###### END OF APPLY CONSTRAINT 
    ########################################################
  }
  
 
  ################################
  ################################
  # RANDOM VARIABLES
  ### Coeficcient matrix
  ################################
  AIC = as.vector((-2 * logL ) + ( 2 * dim(xm)[2]))
  BIC = as.vector((-2 * logL ) + ( log(length(y)) * dim(xm)[2]))
  
  #zvar <- which(unlist(lapply(ZETA, function(x){names(x)[1]})) == "Z")
  
  #####################
  
  varo <- as.list(var.com2)  # variance components for random, no error
  
  Gspo <- lapply(as.list(c(1:length(ZETA))),function(x,K,v){
    oo=K[[x]]*as.numeric((v[[x]])); return(oo)
  }, K=Gs, v=varo) ##
  Gsp <- as(do.call("adiag1", Gspo),Class="sparseMatrix") # in diagonal
  
  Rsp <- R[[1]] * as.numeric(var.com[1 + length(zvar)]) # change
  if(length(R)>1){
    for(tu in 2:length(R)){
      Rsp <- Rsp + as(R[[tu]] *as.numeric(var.com[tu + length(zvar)]),Class="sparseMatrix")
    }
  }
  varo=NULL
  Gspo=NULL
  #########################################################
  # Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467)
  # R-  --  [R-Z[Z'R-Z+G-]-Z'R-]  #-- means minus
  if(sherman){
    Rinv=solve(Rsp,sparse=TRUE, tol = 1e-19)
    Ginv=solve(Gsp,sparse=TRUE, tol = 1e-19)
    ZRZG= solve( tZsp%*%Rinv%*%Zsp + Ginv ,sparse=TRUE, tol = 1e-19)
    vm <- Zsp%*%(Gsp%*%tZsp) + Rsp # V=ZGZ+R
    Vinv = Rinv - ( Rinv%*%Zsp%*%ZRZG%*%t(Zsp)%*%Rinv)
    #########################################################
  }else{
    vm <- Zsp%*%crossprod(Gsp,tZsp) + Rsp # ZGZ+R
    Vinv <- solve(vm,sparse=TRUE, tol = 1e-19)
  }
  #################
  #Vinv2 <- Vinv
  ## Fixed effects
  # B=(X'X V- XVy)-
  xvx <- crossprod(xm, Vinv %*% xm)
  xvxi <- solve(xvx) # variance of fixed effects
  beta <- xvxi %*% crossprod(xm, Vinv %*% y) # (XVX)-XV-y
  #################
  Var.u <- vector(mode="list", length = length(zvar))
  PEV.u <- Var.u
  u <- Var.u
  
  
  pm=Vinv-Vinv%*%xm%*%(xvxi%*%txm%*%(Vinv))
  ZETA3 <- lapply(ZETA, function(x){y=list(Z=as(x[[1]],Class="sparseMatrix"),K=as(x[[2]],Class="sparseMatrix"))})
  
  ee <-  (y - (xm %*% beta))
  
  if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
    for (k in zvar) { # u = KZ'V- (y- XB)
      Uxi <- solve(t(Us[[k]]))
      u[[k]] <- Uxi %*% ( ( (ZETA3[[k]][[2]]*as.numeric(var.com2[k,1])) %*% t(ZETA3[[k]][[1]]) %*% Vinv %*% ee ))
      # (sigma2.u^2) *  Usi  %*% tcrossprod( crossprod(Z%*%K, P)  %*%  (Z%*%K),Usi  )
      Var.u[[k]] <- (as.numeric(var.com2[k,1])^2) *  Uxi %*% tcrossprod( crossprod(ZETA3[[k]][[1]]%*%ZETA3[[k]][[2]], pm)  %*%  (ZETA3[[k]][[1]]%*%ZETA3[[k]][[2]]), Uxi   ) # sigma^4 ZKP ZK
      #Var.u[[k]] <- (as.numeric(var.com2[k,1])^2) *  ( crossprod(ZETA3[[k]][[1]]%*%ZETA3[[k]][[2]], pm)  %*%  (ZETA3[[k]][[1]]%*%ZETA3[[k]][[2]])   ) # sigma^4 ZKP ZK
      PEV.u[[k]] <- as.numeric(var.com2[k,1]) * ZETA3[[k]][[2]] - Var.u[[k]]  # standard errors (SE) for each individual
    }
  }else{
    for (k in zvar) { # u = GZ'V- (y- XB)  
      # where GZ'V- is the proportion of the variance, i.e. %, and then multiplied by the residual or random part leftover after removing fixed effects
      u[[k]] <- ( (ZETA3[[k]][[2]]*as.numeric(var.com2[k,1])) %*% t(ZETA3[[k]][[1]]) %*% Vinv %*% ee )
      Var.u[[k]] <- (as.numeric(var.com2[k,1])^2) *  ( crossprod(ZETA3[[k]][[1]]%*%ZETA3[[k]][[2]], pm)  %*%  (ZETA3[[k]][[1]]%*%ZETA3[[k]][[2]])   ) # sigma^4 ZKP ZK
      PEV.u[[k]] <- as.numeric(var.com2[k,1]) * ZETA3[[k]][[2]] - Var.u[[k]]  # standard errors (SE) for each individual
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
  ####################
  ###################
  ###################
  if(!is.null(names(ZETA))){
    names(u) <- names(ZETA)
    names(Var.u) <- names(ZETA)
    names(PEV.u) <- names(ZETA)
  }
  logL <- as.vector(logL)
  ###################
#   if(!abnormal & !abnormalVE){
#     jkl <- c(23,18,9,20,20,5,14, NA,2,25,NA,7,9,15,22,1,14,14,25,NA,3,15,22,1,18,18,21,2,9,1,19)
#     oh.yeah <- paste(letters[jkl],collapse = "")
#     cat("\nResidual variance (Ve) was pushing to be negative. Be careful, your model might be overspecified.\nPlease check and delete random effects if neccesary.\n")
#   }
  ####
  if(EIGEND){
    Usi <- solve(t(Usp))
    ee <- Usi %*% ee # adjusted residuals
    residuals3 <- (Usi%*%y) - fitted.y[good]
    Vinv <- Usi%*%tcrossprod(Vinv, Usi) # adjusted V inverse
  }
  
  
  out1 <- as.matrix(var.com2, ncol=1); colnames(out1) <- "component" # variance components
  
  rownames(out1) <- c(paste("Var(",varosss,")",sep=""))
  
  
  ### let user knows the constraints
  out1 <- as.data.frame(out1)
  out1[,"constraint"] <- "Positive"
  if(length(boundary)>0){
    out1[boundary,"constraint"] <- "Boundary"
    out1[boundary,1] <- 0
  }
  res <- list(var.comp=out1, V.inv = Vinv, u.hat=u, Var.u.hat=Var.u, 
              PEV.u.hat=PEV.u, beta.hat=beta, Var.beta.hat=xvxi, 
              LL=logL, AIC=AIC, BIC=BIC, X=xm, fitted.y=fitted.y, 
              fitted.u=fitted.u, residuals=ee, cond.residuals=residuals3,
              fitted.y.good=fitted.y.good, Z=Zsp, K=Ksp, fish.inv=aimi,
              convergence=convergence, sigma.scaled=var.com)
  
  layout(matrix(1,1,1))
  #print(big.peaks.col(logL2.stored, -100000))
  return(res)
} # end of else statment at the very beggining to decide if run or not




