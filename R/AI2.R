AI2 <- function(y, X=NULL, ZETA=NULL, R=NULL, draw=TRUE, REML=TRUE, silent=FALSE, init=NULL, iters=50, sherman=FALSE){
  y.or <- y
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  ###########################
  ## y is a vector for the response variable
  ## X is an incidence matrix for fixed effects
  ## Z is a list of lists for each random effect
  # the list of Z can or cannot include the covariance matrix for such random effect
  # if provided must be provided as Z=list(list(Z=Z1,K=K1),list(Z=Z2,K=K2), etc) 
  ############################
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
    if(is.null(R)){R <- diag(length(y))} # dim(x[[1]])[2]
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
    ## to use later for fitted values
    x.or <- as.matrix(xm)
    zeta.or <- ZETA
    zeta.or  <- lapply(zeta.or , function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    ##
    ZETA2 <- ZETA; y2 <- y ; good <- which(!is.na(y)) # make a copy and identify no missing data
    ZETA <- lapply(ZETA2, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
    y <- y[good]
    ZETA <- lapply(ZETA, function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    xm <- as.matrix(xm[good,])
    R <- R[good,good]
    #######################
    ## Initial values
    var.y <- var(y, na.rm=TRUE)
    yv <- scale(y)
    nvarcom <- length(ZETA) + 1
    base.var <- var(yv, na.rm = TRUE)/nvarcom
    if(is.null(init)){var.com <- rep(base.var, nvarcom)}else{var.com<-init/var.y}
    # at the end error variance
    zvar <- which(unlist(lapply(ZETA, function(x){names(x)[1]})) == "Z")
    tn = length(yv)
    logL2=-10000000 # initial log-likelihood
    logL2.stored <- round(logL2,0)
    conv=0 # convergence
    wi=0 # counter
    record<- matrix(var.com*var.y, ncol=1)
    ###########################
    #### BECOME SPARSE
    #R=as(R,Class="sparseMatrix")
    Zs <- lapply(ZETA, function(x){x[[1]]})
    Gs <- lapply(ZETA, function(x){x[[2]]})
    Zsp <- as(do.call("cbind", Zs),Class="sparseMatrix") # column bind Z=[Z1 Z2 Z3]
    fail=FALSE
    # use next for first derivatives in algorithm
    ZETA2 <- lapply(ZETA, function(x){y=list(Z=as(x[[1]],Class="sparseMatrix"),K=as(x[[2]],Class="sparseMatrix"))})
    om <- list()
    for (k in zvar) {
      om[[k]] <- tcrossprod(ZETA2[[k]][[1]], ZETA2[[k]][[1]] %*% (ZETA2[[k]][[2]]) ) 
    }
    om[[length(om)+1]] <- as(diag(length(y)),Class="sparseMatrix") 
    
    ##################
    ## initialize the progress bar
    if(!silent){
      count <- 0
      tot <- 15
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    #####################
    while (conv==0) { # ==================== START =========================
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
      Rsp <- as(R*as.numeric(var.com[length(var.com)]),Class="sparseMatrix") # R matrix
      #########################################################
      # Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467)
      # R-  --  [R-Z[Z'R-Z+G-]-Z'R-]  #-- means minus
      if(sherman){
        Rinv=solve(Rsp)
        Ginv=solve(Gsp)
        ZRZG= solve( as(t(Zsp)%*%Rinv%*%Zsp + Ginv, Class="sparseMatrix")  )
        vm <- Zsp%*%Gsp%*%t(Zsp) + Rsp # V=ZGZ+R
        vmi = Rinv - ( Rinv%*%Zsp%*%ZRZG%*%t(Zsp)%*%Rinv)
        #########################################################
      }else{
        vm <- Zsp%*%Gsp%*%t(Zsp) + Rsp # V=ZGZ+R
        vmi <- solve(vm) # inverse of V
      }
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # P matrix (page 290 in AJHG 96:283-294), projection matrix
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      xvx=t(xm)%*%vmi%*%xm # X'V-X
      xvxi=solve(xvx) # (X'V-X)-
      s1=vmi%*%xm # in steps to make computations faster
      s2=xvxi%*%t(xm)%*%vmi 
      pm=vmi-s1%*%s2
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # log Likelihood (page 290 in AJHG 96:283-294)
      ### ----------------------------------------------------------------- ###
      ### ----------------------------------------------------------------- ###
      # covariance matrices need to be full rank or negative values will indicate matrices that are not positive semidefinite
      if(REML==TRUE){ # ======= ********* WHEN REML ************ ==================
        ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
        if(is.infinite(ddv)){ # ======= IF DETERMINANT IS INFINITE(REML) =======
          stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.")
        }else{ # ======= IF EVERYTHING GOES WELL =========
          logL=as.numeric(-0.5*((ddv)+log(det(xvx))+t(yv)%*%pm%*%yv)) # log likelihood, problem
        }
      }else{ # ======= ********* WHEN ML ************ ==================
        ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
        if(is.infinite(ddv)){ # ======= IF DETERMINANT IS INFINITE (ML) =======
          stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.")
        }else{ # ======= IF EVERYTHING GOES WELL(ML) =========
          logL=as.numeric(-0.5*((ddv) + t(yv)%*%pm%*%yv )) # log likelihood, problem
        }
      }
      #############
      ############# WAS abs(logL-logL2)<0.001     
      if (abs(logL-logL2)<0.001 | wi == iters ){ ## CONVERGENCE, not sure if should be absolute value or not
        conv=1
        if(!silent){
          setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
        }
        if(wi == iters){ # LAST RESOURCE CHANGE TO EM ALGORITHM
          cat("\nMaximum number of iterations reached with no convergence. Changing to EM algorithm\n")
          LRes <- EM(y=y, X=X, ETA=ZETA, R=R, iters = iters, REML=REML, draw=draw, silent=silent)
          #var.com2 <- LRes$var.comp
          logL <- LRes$LL
          fail <- TRUE
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
        aim=matrix(0,nvarcom,nvarcom) # average information matrix
        py=pm%*%yv # Py
        
        for(i in 1:nvarcom){
          for(j in 1:i){
            aim[i,j]=as.numeric(0.5*(t(yv)%*%om[[i]]%*%pm%*%om[[j]]%*%pm%*%py)) # .5 * y' Hi %*% P %*% Hj %*% P %*% P %*% y
            # where Hi=ZKZ, here om[[i]]
            if(i != j){
              aim[j,i] <- aim[i,j]
            }
          }
        }
        aimi=solve(aim) # solve?? (AI)- # Average Information Matrix Inverse (aimi)
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        # D = matrix of first derivatives 
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
          tr1=sum(diag(prm))
          # fill the first derivative matrix
          dldv[k,1]=-0.5*tr1+0.5*as.numeric(t(yv)%*%prm%*%py) # -.5 * [ tr(PHi)+0.5 * y' %*% PHi %*% Py ]        ## fill 1st derivative matrix 1,1
        }
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        # AI update (page 290 in AJHG 96:283-294 or Eq.7 in GSE 38: 25-43)
        up=aimi%*%dldv # updated = (AI)- D
        ### ----------------------------------------------------------------- ###
        ### ----------------------------------------------------------------- ###
        var.com <- as.matrix(var.com) + up
        # CONSTRAINTS: variance components cannot be zero or exceed total variance
        fail <- which(var.com <= 0)
        if(length(fail) > 0){
          var.com[fail] <- 0#base.var/sample(2:4,1)#1e-08#check that now if var.com < 1e-17 is considered zero
        }
        extreme <- which(var.com > 1)
        if(length(extreme) > 0){
          var.com[extreme] <- .002#1e-06#base.var
        }
        record <- cbind(record, var.com * as.numeric(var.y))
        ###########
        lege2 <- list()
        for(k in 1:length(var.com)){
          if(k == length(var.com)){
            lege2[[k]] <- paste("Var(e):")
          }else{
            lege2[[k]] <- paste("Var(u",k,"):",sep="")
          }
        }
        ###########
        # PLOT
        if(draw){
          ylim <- max(unlist(record), na.rm=TRUE)
          my.palette <- brewer.pal(7,"Accent")
          layout(matrix(1:2,2,1))
          plot(logL2.stored[-1],type="l", main="logLikelihood", col=my.palette[7],lwd=3, las=2, xaxt="n", ylab="logLikelihood value", xlab="Iterations processed", cex.axis=0.5) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          plot(record[1,],ylim=c(0,ylim),type="l", las=2, xaxt="n",main="Average Information algorithm results", col=my.palette[1],lwd=3, ylab="Value of the variance component", xlab="Iterations processed", cex.axis=0.5) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          for(t in 1:(dim(record)[1])){
            lines(record[t,],col=my.palette[t],lwd=3)
          } 
          
          ww <- dim(record)[1]
          lege <- list()
          #lege2 <- list()
          for(k in 1:length(var.com)){
            if(k == length(var.com)){
              lege[[k]] <- paste("Var(e):",round(record[k,wi+1],4), sep="")
              #lege2[[k]] <- paste("Var(e):")
            }else{
              lege[[k]] <- paste("Var(u",k,"):",round(record[k,wi+1],4), sep="")
              #lege2[[k]] <- paste("Var(u",k,"):",sep="")
            }
          }
          legend("topleft",bty="n", cex=0.7, col=my.palette, lty=1, lwd=3, legend=unlist(lege))
        }
        ###########
        fail=FALSE
        ###########
      }
      
    }# final variance components here:
    if(fail){
      var.com2 <- as.matrix(LRes$var.comp)
    }else{
      var.com2 <- as.matrix(record[,dim(record)[2]])
    }
  }
  return(var.com2)
}


