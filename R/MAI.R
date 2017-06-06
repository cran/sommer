MAI <- function(Y, X=NULL, ZETA=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=20, init=NULL, tol=1e-3, che=TRUE, EIGEND=FALSE, forced=NULL, IMP=FALSE){
  
  Yh <- Y
  Xh <- X
  Zh <- ZETA
  
  if(tol == 1e-22){
    MARIA <- function(y, X=NULL, ZETA=NULL, R=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=50, constraint=TRUE, init=NULL, sherman=FALSE, che=TRUE, EIGEND=FALSE, Fishers=FALSE, gss=TRUE, forced=NULL){
      
      if(EIGEND){
        DISO <- dim(ZETA[[1]]$Z)
        if(DISO[1] != DISO[2]){
          stop("EIGEN DECOMPOSITION EIGEND ONLY WORKS FOR SQUARE PROBLEMS 
               'Z' MATRIX IS NOT SQUARE",call.=FALSE)
          
        }
      }
      
      y.or <- y
      x.or <- X
      pushess <- rep(0,length(ZETA)+1)
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
        if(is.null(R)){R <- diag(length(y))} # dim(x[[1]])[2]
        
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
        xm <- as.matrix(xm[good,])
        txm <- t(xm)
        R <- R[good,good]
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
        om[[length(om)+1]] <- as(diag(length(y)),Class="sparseMatrix") 
        #######################
        ## Initial values
        if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])){
          y <- as.vector(t(Usp) %*% as.matrix(y,ncol=1))
          xm <- t(Usp) %*% xm
          txm <- t(xm)
        }
        var.y <- var(y, na.rm=TRUE)
        yv <- scale(y)
        nvarcom <- length(ZETA) + 1
        base.var <- var(yv, na.rm = TRUE)#/nvarcom
        ### EIGEND requires very small initial var.comp compared to regular models
        if(length(ZETA)==1 & EIGEND==TRUE & (dim(ZETA[[1]][[1]])[1] == dim(ZETA[[1]][[2]])[2])){
          if(is.null(init)){var.com <- c(rep(.0001, nvarcom))}else{var.com <- init/var.y} # at the end error variance
          #print(var.com)
        }else{ # before instead of base.var was .01
          if(is.null(init)){var.com <- c(rep(base.var, nvarcom))}else{var.com <- init/var.y} # at the end error variance
        }
        weird=FALSE
        tn = length(yv)
        logL2=-10000000 # initial log-likelihood
        logL2.stored <- round(logL2,0)
        conv=0 # convergence
        wi=0 # counter
        record<- matrix(var.com*var.y, ncol=1)
        taper <- rep(0.9, iters) # weighting parameter for updates
        taper[1:2] <- c(0.5, 0.7)
        #var.com <- rep(var(yy2, na.rm = TRUE)/nvarcom, nvarcom) ### NEW at the end error variance
        ##################
        if(is.null(names(ZETA))){
          varosss <- c(paste("u.",df.ord, sep=""))
        }else{
          varosss <- c(names(ZETA))
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
        #gss <- TRUE
        #sumdiags <- unlist(lapply(ZETA,function(x){is.diagonal.matrix(x[[2]])}))
        #if(length(which(sumdiags)) == length(ZETA)){gss <- FALSE}
        ######################
        #var.com <- NR2(y=y, X=xm, ZETA=ZETA, R=R, REML=REML, draw=FALSE, silent=TRUE, iters=2)
        
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
        
        ups <- numeric()
        if(is.null(forced)){ # &&&&&&&&&&&&&&&& IF NOT FORCED &&&&&&&&&&&&&&&&&
          
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
            Rsp <- as(R*as.numeric(var.com[length(var.com)]),Class="sparseMatrix") # R matrix
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
              vmi <- solve(vm, sparse=TRUE, tol = 1e-19) # inverse of V
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
            #print(pm[1:5,1:5])
            #print(var.com)
            #vmi.xm <- vmi%*%xm
            #xvx=txm%*%vmi.xm # X'V-X
            #xvxi=solve(xvx, sparse=TRUE, tol = 1e-19) # (X'V-X)-
            #s1=vmi%*%xm # in steps to make computations faster
            #s2=xvxi%*%txm%*%vmi 
            #pm=vmi-crossprod(t(vmi.xm),s2) #vmi-s1%*%s2#
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
            if(wi > 10 & gss == TRUE){
              
              joke <- round(dim(ups)[2]*.2)
              todo <- (dim(ups)[2]-joke):dim(ups)[2]
              pushes <- apply(ups[,todo],1,function(fg){length(which(fg < 0))/length(fg)}) # how many times were pushed to be negative
              kkkkk <- length(which(pushes > .4))
              
              if( ( zig.zag(logL2.stored[(length(logL2.stored)-4):length(logL2.stored)]) == 1) & (kkkkk == 0)  ){
                wi=iters
                if(!silent){
                  setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
                }
                cat("\nA weird likelihood behavior has been encountered. \nBe careful with variance components returned.\nSystem has singularities, ML estimators returned.\nPlease check results trying the NR or EM algorithms")
                draw=TRUE
              }else if(( zig.zag(logL2.stored[(length(logL2.stored)-4):length(logL2.stored)]) == 1) & (kkkkk != 0)){
                #print(pushes)
                wi=iters
                
              }
              
            }
            #################################
            
            #############
            ############# WAS abs(logL-logL2)<0.001     
            if (abs(logL-logL2)<0.001 | wi == iters ){ ## CONVERGENCE, not sure if should be absolute value or not
              conv=1
              if(!silent){
                setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
              }
              if(wi == iters){ # LAST RESOURCE CHANGE TO EM ALGORITHM
                ## zig.zag function
                #picos <- (big.peaks.col(logL2.stored, -10000000))
                #picos <- list(pos=picos$pos[which(picos$pos > 3)], hei=picos$hei[which(picos$pos > 3)])
                #MLE <- picos$pos[1]
                MLE <- which(logL2.stored == max(logL2.stored))[1]
                ##
                last20 <- round(length(logL2.stored)*.2)
                best.try <- record[,(dim(record)[2]-last20):(dim(record)[2])]
                
                logL2.stored[which(logL2.stored == max(logL2.stored))]
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
              #print(dldv)
              ### ----------------------------------------------------------------- ###
              ### ----------------------------------------------------------------- ###
              # AI update (page 290 in AJHG 96:283-294 or Eq.7 in GSE 38: 25-43)
              
              up=aimi%*%dldv
              #print(up)
              ### ----------------------------------------------------------------- ###
              ### ----------------------------------------------------------------- ###
              
              failup <- which(up > 0.2)
              if(length(failup) > 0){
                up[failup,] <- up[failup,]/5#mean(up[-failup,],na.rm=TRUE)#var.com[fail]/2#1e-3#base.var/sample(2:4,1)#1e-08#check that now if var.com < 1e-17 is considered zero
              }
              #print(up)
              
              #var.com <- as.matrix(var.com) + as.matrix(up)
              var.com <- as.matrix(var.com) + (taper[wi]*as.matrix(up))# taper is a regress trick
              #var.como <- var.com;
              #colnames(var.como)<-"A"
              #print(var.como)
              #print(up)
              ### if likelihood starts to take a weird shape and is an mmer2 model
              if(((abs(logL2.stored[length(logL2.stored)]) - logL2.stored[length(logL2.stored)-1]) > 0) & (wi > 15) & (gss == FALSE)){
                non.zero <- which(as.vector(up) < -.25)
                if(length(non.zero) > 0){var.com[non.zero] <- base.var}
                
                LRes <- EM2(y = y, X = X, ETA = ZETA, R = R, 
                            iters = 5, REML = REML, draw = FALSE, 
                            silent = silent, init=var.com*var.y)
                var.com <- as.vector(LRes)/var.y
                weird=TRUE
              }
              
              ups <- cbind(ups,var.com)
              #print(var.com)
              # CONSTRAINTS: variance components cannot be zero or exceed total variance
              # if after 5 iterations keeps pushing to negative values make it zero
              fail <- which(var.com <= 0)
              if(length(fail) > 0){
                if(nmbb){
                  var.com[fail] <- .05
                }else{
                  var.com[fail] <- 1e-3
                }
                #var.com[fail]/2#1e-3#base.var/sample(2:4,1)#1e-08#check that now if var.com < 1e-17 is considered zero
              }
              extreme <- which(as.vector(unlist(var.com)) > 1.5)
              if(length(extreme) > 0 & wi > 1){ # just do it after the second iteration
                #stn <- max(var.com[-extreme])+.05
                if(nmbb){
                  var.com[extreme] <- .05
                }else{
                  var.com[extreme] <- 1e-3
                }
                #record[extreme,(wi-1)]/var.y#.01#var.com[extreme]/2# 2e-3#1e-06#base.var
              }
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
                #lege2 <- list()
                for(k in 1:length(var.com)){
                  if(k == length(var.com)){
                    lege[[k]] <- paste("Var(e):",round(record[k,wi+1],4), sep="")
                    #lege2[[k]] <- paste("Var(e):")
                  }else{
                    lege[[k]] <- paste("Var(",varosss[k],"):",round(record[k,wi+1],4), sep="")
                    #lege2[[k]] <- paste("Var(u",k,"):",sep="")
                  }
                }
                legend("topleft",bty="n", cex=0.7, col=my.palette, lty=1, lwd=3, legend=unlist(lege))
              }
              
              ##$$$$$$$$$$$$$$$$
              ##### June 1 2016
              ### sometimes the variance component gets pushed to zero
              ### but never achieves the zig zag likelihood or pushed enough to be considered zero
              ## although is actually zero
              abnormal <- FALSE
              abnormalVE <- TRUE
              if(wi > 10){
                joke <- round(dim(ups)[2]*.2)
                todo <- (dim(ups)[2]-joke):dim(ups)[2]
                pushes <- apply(ups[,todo],1,function(fg){length(which(fg < 0))/length(fg)}) # how many times were pushed to be negative
                pushess <- rbind(pushess,pushes)
                crush <- apply(as.matrix(pushess),2,sum)
                badbad<- which(crush > 1)
                #print("YES")
                #print(badbad);print(abnormal)
                if((length(badbad)>0) ){ # if there's a negative component
                  wi=iters-1
                  # error IS NOT in the bad
                  if(length(crush) %!in% badbad){ # if a random effect other than error variance is negative
                    abnormal<- TRUE 
                  }else{ # if error variance is negative IS IN THE BAD
                    if(constraint){
                      cat("\nModel is overspecified. A model with no error variance (or negative) is not reliable.\nReturning maximum likelihood estimators. Optionally, you could remove some random \neffects or try a different algorithm.\n")
                      constraint=FALSE
                      #abnormalVE <- FALSE
                      #stop()
                    }else{
                      #cat("\nError variance (Ve) was pushing to be negative. You might be overspecifying your model.")
                      abnormalVE <- FALSE # means the opposite that abnormal VE was found but is faster to check for TRUE IF statements
                    }
                  }
                }
              }
              #print(abnormal)
              ##$$$$$$$$$$$$$$$$
              
              ###########
              fail=FALSE
              ###########
            }
            
          }# =====================================  END =======================================
          #####################################################################################
          ####### FOR LOOP FOR ITERATIONS FINISHED ############################################
          if(fail){
            var.com2 <- as.matrix(record[,MLE])
          }else{
            var.com2 <- as.matrix(record[,dim(record)[2]])
          }
          #########################################################
          # perform boosting if there was a var comp close to zero
          # and AI is able to converge
          
          #just in the last iterations we check if still trying to be zero
          joke <- round(dim(ups)[2]*.2)
          todo <- (dim(ups)[2]-joke):dim(ups)[2]
          
          pushes <- apply(ups[,todo],1,function(fg){length(which(fg < 0))/length(fg)}) # how many times were pushed to be negative
          nnn <- length(which(pushes >= .4));
          if(length(nnn) > 0){fail=FALSE}
          #print(pushes);print(ups); print(fail); print(nnn);print(abnormal);print(abnormalVE)
          zero <- which(pushes >= 0.4 )
          
          if(length(pushes) %in% zero){
            abnormal=FALSE 
            abnormalVE=FALSE
          }
          
          RE <- length(ZETA) # how many random effects
          
          #print(RE)
          
          #print(constraint);print(EIGEND);print(abnormal)
          #print(pushes)
          ###### APPLY CONSTRAINT 
          ###### APPLY CONSTRAINT 
          if( (RE > 1)){ #only apply if there is more than one variance component
            
            if(length(ZETA) != nnn){
              
              if(((nnn >= 1 & nnn <= 4) & (nnn < (dim(var.com2)[1])-1) & (fail == FALSE) & (constraint == TRUE) & (abnormalVE)) | abnormal ){
                
                if(!silent){
                  cat("\nOne or more variance components pushing to be zero. Boundary constraint applied.\n")
                }
                zero <- which(pushes >= 0.4 ) # find zero var.comp and remove
                if(abnormal){
                  zero <- badbad
                }
                nonzero <- (1:dim(var.com2)[1])[-zero]
                #print(zero);print(nonzero);print(abnormal)
                ## estimate accurately the good variance components
                allfailed <- zero==1:(length(ZETA))
                if(length(which(allfailed)) == length(ZETA)){ # if all var.comp are zero start all over again with a different algorithm
                  ##sometimes imputing the missing data when GWAS leads to very few information causing all variance
                  ## components to be zero and AI cannot do anything, we change to NR
                  boost <- NR22(y=y.or, X=x.or, ZETA=zeta.or, R=NULL, REML=REML, draw=draw, silent=silent, iters=20, sherman=sherman)
                  var.com2[,1] <- boost
                }else{
                  boost <- ai2help(y=y.or, X=x.or, ZETA=zeta.or[-zero], R=NULL, REML=REML, draw=draw, silent=silent, iters=20, init=as.vector(var.com2)[-zero], sherman=sherman)
                  var.com2[nonzero,] <- boost
                  var.com2[zero,] <- 0
                }
                
                #var.com2[zero,] <- 5e-5
                ## force found values and get the EM closest to those values
                #boost2 <- AI3(y=y.or, X=x.or, ZETA=zeta.or, R=NULL, REML=REML, draw=draw, forced=nonzero, silent=silent, iters=15, init=as.vector(var.com2/var.y),sherman=sherman)
                #boost2 <- EM2(y=y.or, X=x.or, ETA=zeta.or, R=NULL, REML=REML, draw=draw, forced=nonzero, silent=FALSE, iters=10, init=as.vector(var.com2))
                #var.com2[zero,] <- boost2[zero,]
              }
              
            }else{
              zero <- which(pushes >= 0.4 ) 
              if(abnormal){
                zero <- badbad
              }
              var.com2[zero,] <- 0  
            }
            
          }else{ #only one var.comp but there was a zero
            if(((nnn >= 1 & nnn <= 4) & (nnn < (dim(var.com2)[1])-1) & (fail == FALSE) & (constraint == TRUE) & (abnormalVE)) | abnormal ){
              zero <- which(pushes >= 0.4 ) 
              if(abnormal){
                zero <- badbad
              }
              var.com2[zero,] <- 0
            }
          }
          ###### ENF OF APPLY CONSTRAINT 
          ###### END OF APPLY CONSTRAINT 
          ########################################################
        }else{ ## &&&&&&&&&&&&  "IF FORCED" &&&&&&&&&&&&&&
          abnormal=FALSE
          abnormalVE=TRUE
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
        ###########################
      } # end of else statment at the very beggining to decide if run or not
      
      
      if(weird){
        cat("\nWeird likelihood behavior found, Additional EM steps were taken to get better initial values\n")
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
      Rsp <- as(R*as.numeric(var.com2[length(var.com2)]),Class="sparseMatrix")
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
      #xvx=t(xm)%*%Vinv%*%xm # X'V-X
      #xvxi=solve(xvx) # (X'V-X)-
      #Vinv.xm <- Vinv%*%xm
      #xvx=txm%*%Vinv.xm # X'V-X
      #xvxi=solve(xvx, sparse=TRUE, tol = 1e-19) # (X'V-X)-
      #s1=Vinv%*%xm # in steps to make computations faster
      #s2=xvxi%*%txm%*%Vinv 
      #pm=Vinv-crossprod(t(Vinv.xm),s2) #Vinv-s1%*%s2#
      
      pm=Vinv-Vinv%*%xm%*%(xvxi%*%txm%*%(Vinv))
      ZETA3 <- lapply(ZETA, function(x){y=list(Z=as(x[[1]],Class="sparseMatrix"),K=as(x[[2]],Class="sparseMatrix"))})
      #for(h in zvar){
      #  
      #  if(EIGEND){
      #    Var.u[[h]] <- (as.numeric(var.com2[h,1])^2) *  ( crossprod(ZETA3[[h]][[1]]%*%ZETA3[[h]][[2]], pm)  %*%  (ZETA3[[h]][[1]]%*%ZETA3[[h]][[2]])   ) # sigma^4 ZKP ZK
      #    PEV.u[[h]] <- as.numeric(var.com2[h,1]) * ZETA3[[h]][[2]] - Var.u[[h]]  # standard errors (SE) for each individual
      #  }else{
      #    
      #  }
      #  
      #} #(sigma2.u^2) *  ( crossprod(Z%*%K, P)  %*%  (Z%*%K)   )
      #################
      ## Random effects
      #u <- list() # we put it up
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
      #print(abnormal); print(abnormalVE)
      if(!abnormal & !abnormalVE){
        jkl <- c(23,18,9,20,20,5,14, NA,2,25,NA,7,9,15,22,1,14,14,25,NA,3,15,22,1,18,18,21,2,9,1,19)
        oh.yeah <- paste(letters[jkl],collapse = "")
        cat("\nResidual variance (Ve) was pushing to be negative. Be careful, model might be overspecified.\n")
      }
      ####
      if(EIGEND){
        Usi <- solve(t(Usp))
        ee <- Usi %*% ee # adjusted residuals
        residuals3 <- (Usi%*%y) - fitted.y[good]
        Vinv <- Usi%*%tcrossprod(Vinv, Usi) # adjusted V inverse
      }
      
      
      out1 <- as.matrix(var.com2, ncol=1); colnames(out1) <- "Variance Components" # variance components
      rownames(out1) <- c(paste("Var(",varosss,")",sep=""), "Var(Residual)")
      res <- list(var.comp=out1, V.inv = Vinv, u.hat=u, Var.u.hat=Var.u, 
                  PEV.u.hat=PEV.u, beta.hat=beta, Var.beta.hat=xvxi, 
                  LL=logL, AIC=AIC, BIC=BIC, X=xm, fitted.y=fitted.y, 
                  fitted.u=fitted.u, residuals=ee, cond.residuals=residuals3,
                  fitted.y.good=fitted.y.good, Z=Zsp, K=Ksp, fish.inv=aimi)
      
      layout(matrix(1,1,1))
      #print(big.peaks.col(logL2.stored, -100000))
      return(res)
      }
    
    
    
    }
  ## identify missing data across variables
  if(IMP){
    ytouse <- 1:nrow(Y)
  }else{
    nona <- unlist(apply(Y,2,function(x){which(!is.na(x))}))
    names(nona) <- NULL
    nonat <- table(nona)
    ytouse <- as.numeric(names(nonat)[which(nonat == dim(Y)[2])])
    ## create provisional data sets for estimating variance components
    Y <- Y[ytouse,]
    if(!is.null(X)){
      X <- X[ytouse,] 
    }else{Xh <- X}
    ZETA <- lapply(ZETA, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=ytouse)
  }
  
  if(EIGEND){
    DISO <- dim(ZETA[[1]]$Z)
    if(DISO[1] != DISO[2]){
      stop("EIGEN DECOMPOSITION EIGEND ONLY WORKS FOR SQUARE PROBLEMS 
           'Z' MATRIX IS NOT SQUARE",call.=FALSE)
      
    }
    cat("EIGEND feature activated. Eigen decomposition of K will be performed\n")
  }
  
  namesY <- colnames(Y)
  if(is.null(namesY)){
    namesY <- paste("T",1:dim(as.matrix(Y))[2],sep="")
  }
  Y <- as.data.frame(Y)
  ###########################
  ### define useful functions
  best.layout <-function(x){
    x1 <- merge(1:x, 1:x)
    x2 <- abs(apply(x1, 1, function(x, des) {
      x <- unlist(x)
      y <- (x[[1]] * x[[2]]) - des
      return(y)
    }, des = x))
    vw <- x1[which(x2 < 2), ]
    x3 <- abs(apply(vw, 1, function(x) {
      y <- abs(x[1] - x[2])
      return(y)
    }))
    x4 <- vw[which(x3 == min(x3))[1], ]
    res <- as.vector(unlist(x4))
    return(res)
  }
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
  ### useful functions defined
  ###########################
  
  ##############################
  ######## CONTROLS ############
  ##############################
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
    #### now check dimensions
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
  
  if(is.null(X) & is.null(ZETA)){
    stop("Empty model",call. = FALSE)
  }else if(is.null(X) & !is.null(ZETA)){
    X <- matrix(rep(1,dim(as.data.frame(Y))[1]))
    if(length(ZETA)==1 & EIGEND==TRUE){
      EIGENS <- lapply(ZETA, function(x){eigen(x[[2]])}) # eigen decomposition of K
      Us <- lapply(EIGENS, function(x){x$vectors}) # extract eigen vectors U
      Usp <- as(do.call("adiag1", Us),Class="sparseMatrix") # U'G as diagonal
      Ds <- lapply(EIGENS, function(x){diag(x$values)}) # extract eigen values D
      Dsp <- as(do.call("adiag1", Ds),Class="sparseMatrix") # U'G as diagonal
      ZETA <- lapply(as.list(1:length(ZETA)),function(x,zz,kk){list(Z=zz[[x]][[1]], K=kk[[x]])}, zz=ZETA, kk=Ds)
      Y <-apply(Y,2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
      Y <- as.matrix((t(Usp) %*% as.matrix(Y)))
      X <- as.matrix(t(Usp) %*% X)
      tX <- t(X)
    }
    #firstcheck <- unlist(lapply(ZETA,function(x){a <- dim(x$Z);if(a[1]>a[2]){STRUC=FALSE}else{STRUC=TRUE};return(STRUC)}))
    #first <- length(which(firstcheck))/length(firstcheck)
    #if(first >=0.5){
    STRUCT=TRUE
    #STRUCT=FALSE
    #}else{STRUCT=FALSE}
    if(EIGEND){
      ZETA <- lapply(ZETA,function(x){x$Z <- as(x$Z,Class="sparseMatrix");x$K <- as(x$K,Class="sparseMatrix"); return(x)}) 
    }else{
      ZETA <- lapply(ZETA,function(x){x$Z <- as(x$Z,Class="sparseMatrix"); return(x)})
    }
    
  }
  if(is.null(names(ZETA))){
    varosss <- c(paste("u.",1:length(ZETA), sep=""),"Residual")
    varosss2 <- c(paste("u.",1:length(ZETA), sep=""))
  }else{
    varosss <- c(names(ZETA),"Residual")
    varosss2 <- c(names(ZETA))
  }
  ##############################
  ######## END CONTROLS#########
  ##############################
  #cat(paste("\n",dim(Y)));cat(dim(X))
  dimos <- dim(as.matrix(Y))
  ts <- dimos[2] #no. of traits
  inds <- dimos[1] #no. of individuals
  nvarcom <- length(ZETA) #no. of var. components
  Y.or <- Y
  Y.or2 <-apply((Y.or),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  Y.or3 <- as.matrix(as.vector(Y.or2)); dim(Y.or3)
  Yt.or3 <- t(Y.or3)
  #cat(dim(Y.or3))
  # impute phenotypes
  Y <-apply((scale(Y)),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  # decompose in a vector
  Y <- as.matrix(as.vector(Y)); dim(Y)
  tY <- t(Y); dim(tY) #get transpose
  X.or <- X
  X <- do.call("adiag1", rep(list(X), ts)); dim(X)
  tX <- t(X); dim(tX)
  qr <- qr(X)
  rankX <- dim(X)[1]-qr$rank
  ## define variance-covariance components, each list is a var.comp
  or.var <- var(Y.or,na.rm=TRUE) #original variance values
  sc.var <- var(scale(Y.or),na.rm=TRUE) #scaled variances
  
  ##initial values
  #for all variance components plus error define initial var.comp values
  if(is.null(init)){
    if(EIGEND){#/(nvarcom+16)
      var.com <- rep(list(sc.var),nvarcom+1)
    }else{#/(nvarcom+6)
      var.com <- rep(list(sc.var/(nvarcom+6)),nvarcom+1)
    }
    
  }else{
    if(is.list(init)){
      ## control of good dimensions
      lapply(init,function(x){if( (dim(x)[1] != dimos[2]) | (dim(x)[2] != dimos[2]) ){
        stop(paste("The customized var-covar values have to be presented in trait x trait matrix.\nIn your case a list with",nvarcom+1,"elements (random effects) storing a",dimos[2],"x",dimos[2],"matrix, not",dim(x)[1],"x",dim(x)[2]),
             call. = FALSE)}})
      var.com <- lapply(init, function(x,y,z){((x*y)/x)/z},y=sc.var,z=nvarcom+2)
      #print("yes")
    }else{
      if( (dim(init)[1] != dimos[2]) | (dim(init)[2] != dimos[2]) ){
        stop(paste("The customized var-covar values have to be presented in trait x trait matrix.\nIn your case a list with",nvarcom+1,"elements (random effects) storing a",dimos[2],"x",dimos[2],"matrix, not",dim(init)[1],"x",dim(init)[2]),
             call. = FALSE)}
      var.com <- rep(list(((init*sc.var)/init)/(nvarcom+2)),nvarcom+1)
    }
  }
  logL2=-10000000
  conv=0
  wi=0
  logLL <- numeric()
  AIsing <- FALSE
  Vsing <- FALSE
  bad2 <- numeric()
  taper <- rep(0.9, iters) # weighting parameter for updates
  taper[1:2] <- c(0.5, 0.7)
  ### ================================================= ###
  ###         Get Vi or dV/ds
  ### ================================================= ###
  
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
  #str(posmats.list.vc)
  ##each element are the derivatives of each var.comp
  deriv.list.vc <- list() 
  for(v in 1:(nvarcom+1)){ ## FOR EACH VAR.COMP
    if(v <= nvarcom){#normal
      deriva <- tcrossprod(ZETA[[v]]$Z %*% ZETA[[v]]$K, (ZETA[[v]]$Z)) #dVi
    }else{#error
      deriva <- diag(dim(ZETA[[1]]$Z)[1])#tcrossprod(ZETA[[1]]$Z %*% diag(dim(ZETA[[1]]$K)[2]), (ZETA[[1]]$Z)) #dVi  #dVi
    }
    #possibles <- posmats.list.vc[[v]]
    if(STRUCT){
      deriv.list.vc[[v]] <- lapply(posmats.list.vc[[v]],function(x,y){as(kronecker(x,y),Class="sparseMatrix")},y=as.matrix(deriva))
    }else{
      deriv.list.vc[[v]] <- lapply(posmats.list.vc[[v]],function(x,y){as(kronecker(x,y),Class="sparseMatrix")},y=as.matrix(deriva))#deriva=ZKZ
    }
    
  } #pos.mats # str(deriv.list.vc) #length(deriv.list.vc)
  ## ================= START ALGORITHM ==================== ##
  ## ================= START ALGORITHM ==================== ##
  ## ================= START ALGORITHM ==================== ##
  ## ================= START ALGORITHM ==================== ##
  ## ================= START ALGORITHM ==================== ##
  if(is.null(forced)){
    ##################
    ## initialize the progress bar
    if(!silent){
      count <- 0
      tot <- 15
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    #####################
    var.coms <- lapply(var.com,function(x){
      if(dim(as.matrix(x))[1]>1){
      base <- upper.tri(x)
    diag(base) <- TRUE
    base2 <- which(base,arr.ind = TRUE)
    base2 <- base2[ order(base2[,1], base2[,2]), ]
    return(as.matrix(x[base2]))}else{return(as.matrix(x))}})
    #print(var.com)
    ups <- do.call("rbind",var.coms)
    ups <- as.matrix(ups*0)
    
    while (conv==0) {
      wi=wi+1
      #####################
      if(!silent){
        count <- count + 1
      }
      #####################
      ### ================================================= ###
      ###          FORM "V" MATRIX 
      ### ================================================= ###
      listGs <- list()
      disk <- lapply(ZETA,function(x){dim(x$K)[2]})
      disk$error <- NA; disk$error <- disk[[1]]
      
      # system.time(asd <-tcrossprod(ZETA[[i]]$Z %*% ZETA[[i]]$K, (ZETA[[i]]$Z))) # same than
      # system.time(asd <- ZETA[[i]]$Z%*%crossprod(ZETA[[i]]$K ,t(ZETA[[i]]$Z)))
      # emmreml solve(kronecker(ZKZt, Vgt) + kronecker(diag(n), Vet) + tolparinv * diag(d * n))
      if(STRUCT){ # kronecker(var.com[[i]],zkz)
        for(i in 1:(nvarcom+1)){ #both ZKZ'
          if(i<=nvarcom){
            zkz <- tcrossprod(ZETA[[i]]$Z %*% ZETA[[i]]$K, (ZETA[[i]]$Z))
            listGs[[i]] <- kronecker(var.com[[i]],as.matrix(zkz))
          }else{
            zkz <- tcrossprod(ZETA[[1]]$Z %*% diag(disk[[i]]), (ZETA[[1]]$Z)) #
            listGs[[i]] <- kronecker(var.com[[i]],as.matrix(zkz))
          }
        }
      }else{ # kronecker(zkz,var.com[[i]])
        for(i in 1:(nvarcom+1)){
          if(i<=nvarcom){ # ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
            zkz <- ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
            listGs[[i]] <- kronecker(as.matrix(zkz),var.com[[i]])
          }else{
            zkz <- diag(dimos[1]) #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
            listGs[[i]] <- kronecker(as.matrix(zkz),var.com[[i]])
          }
        }
      }
      
      V <- matrix(0,dimos[1]*dimos[2],dimos[1]*dimos[2])
      for(l in 1:length(listGs)){
        V <- V + listGs[[l]]
      }; 
      (V[1:5,1:5]);dim(V)
      
      Vmi <- try(solve(as(V, Class="sparseMatrix")), silent = TRUE)
      if (class(Vmi) == "try-error") {
        diag(V) <- diag(V) + (rep(1e-06, dim(V)[2]))
        Vsing <- TRUE
        Vmi <- solve(V)
      }
      dim(V) #dimens
      Vmi[1:5,1:5]
      
      ### ================================================= ###
      ###     P matrix (page 290 in AJHG 96:283-294)
      ### ================================================= ###
      xvx <- crossprod(X, Vmi %*% X)
      P <- Vmi - Vmi %*% X %*% solve(xvx, crossprod(X, Vmi))
      Vmi <- NULL #release memory
      ytPy <- tY%*%(P%*%Y)
      
      var.com <- lapply(var.com,function(x){x*(as.numeric(ytPy)/as.numeric(rankX))})
      P <- P * (as.numeric(rankX)/as.numeric(ytPy))
      ### ================================================= ###
      ###     log Likelihood (page 290 in AJHG 96:283-294)
      ### ================================================= ###
      ddv <- determinant(V, logarithm = TRUE)$modulus[[1]]
      V <- NULL #release memory
      logL=as.numeric(-0.5*((ddv)+determinant(xvx, logarithm = TRUE)$modulus[[1]]+ytPy)) # log likelihood, problem
      #cat("iteration #",wi,logL,"\n")
      logLL[wi] <- logL
      
      if (((logL-logL2 < tol) | (wi==iters)) ) {
        #if(wi > 8){
        conv=1
        if(!silent){
          setTxtProgressBar(pb, (tot-1/tot))### keep filling the progress bar
        }
        #}
        
      }else{
        #if(wi > 8){
        if(!silent){
          setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
        }
        #}
      }
      logL2=logL
      #str(deriv.list.vc)
      #print(STRUCT)
      Py <- crossprod(P,Y)
      ### ================================================= ###
      ###         Fill AI matrix (2nd derivatives)
      ### ================================================= ###
      #make a single list
      at <- do.call(list, unlist(deriv.list.vc, recursive=FALSE))#lapply(at,dim)
      (AI <- matrix(NA, length(at), length(at)))
      for(h1 in 1:dim(AI)[2]){
        for(h2 in 1:h1){
          ##AI 
          AI[h1,h2] <- 0.5 * (tY %*% at[[h1]] %*% P %*% at[[h2]] %*% crossprod(P, Py))[1]
          ##NR
          #term <- at[[h1]] %*% P %*% at[[h2]] %*% P
          #AI[h1,h2] <- 0.5 * (sum(diag(term)) - tY %*% term %*% Py)
        }
      }
      AI <- copying2(AI)
      # AI[1:5,1:5];  
      #image(AI)
      
      AIi <- try(solve(AI), silent = TRUE)
      if (class(AIi) == "try-error") {
        diag(AI) <- diag(AI) + (rep(1e-03, dim(AI)[2]))
        AIsing <- TRUE
        AIi <- try(solve(AI), silent = TRUE)
        if (class(AIi) == "try-error") {
          AIi <- ginv(AI)
        }
      } #AIi[1:5,1:5]
      ### ================================================= ###
      ###          MATRIX OF CURRENT VALUES
      ### ================================================= ###
      # get diagonal and upper triangulars for var.comp matrices
      varos <- lapply(var.com, function(x){
        if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
      babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
      return(x[babas])}else{return(as.matrix(x))}})
      # get indices to know who are variances and who are covariances
      varos2 <- lapply(var.com, function(x){
        if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
      babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
      return(babas)}else{return(cbind(1,1))}})
      
      vcv <- do.call("rbind",varos2) # order of random effects
      
      var.comM <- as.matrix(unlist(varos)) # current values of vc
      
      ### ================================================= ###
      ###          MATRIX OF FIRST DERIVATIVES
      ### ================================================= ### 
      dldv <- numeric()
      for(t1 in 1:length(at)){
        prm=P%*%at[[t1]] # PHi# ( P(e) = Vinv - [Vinv X (X'V-X)- X Vinv]  )  %*% I or K, etc
        tr1=sum(diag(prm)) # trace
        dldv[t1] = -0.5*tr1+0.5*as.numeric(tY%*%prm%*%Py)
      }
      dldv <- as.matrix(dldv)
      ### ================================================= ### 
      ### update (page 290 in AJHG 96:283-294 or Eq.7 in GSE 38: 25-43)
      ### ================================================= ### 
      are.var <- which(vcv[,1]==vcv[,2]) 
      are.covar <- which(vcv[,1] != vcv[,2]) 
      
      up=AIi%*%dldv
      
      #### controled update strategy
      if(wi >=1){ #after the second iteration because will help us to find the var.comp==0
        failup <- which(abs(up) > 0.25)
        # ========================
        # variances cannot be negative, 
        # this portion helps to control them 
        # by identifying and excluding them of the process (only 1st round)
        if(wi==1 & STRUCT){
          #bad <- which((up) < -0.2)
          #bad2 <- intersect(are.var,bad)
        }
        # ============================
        if(length(failup) > 0){
          #cat(paste("\n",up[failup,]))
          up[failup,] <- ups[failup,wi]*.7
        }
      }
      #print(dim(up))
      #print(dim(ups))
      ups <- cbind(as.matrix(ups),as.matrix(up))
      #### end controled update strategy
      #print(ups)
      
      # update
      #theta <- var.comM + up # average information  #theta <- var.comM - up # newton-raphson
      theta <- var.comM + (taper[wi]*up)
      #print(theta)
      # extreme values with scaled variance > 1 are illegal
      ext <- which(abs(theta) > 1.1)
      if(length(ext) > 0){
        #bad2 <- intersect(are.var,ext)
        #bad2 <- vcv[ext,]
        #theta[ext,] <- 0.5
      }
      
      ## see the number of variance compnents to extract for each, i.e. add, dom, error
      no.var <- lapply(deriv.list.vc,function(x){length(x)})
      
      for(r in 1:length(no.var)){
        si <- 1:no.var[[r]]
        newmat <- matrix(NA,ts,ts)
        sq <- upper.tri(newmat)
        diag(sq) <- TRUE
        babas2 <- which(sq,arr.ind = TRUE)
        babas2 <- babas2[ order(babas2[,1], babas2[,2]), ]
        newmat[babas2] <- theta[si,1]
        var.com[[r]] <- copying(newmat)
        theta <- matrix(theta[-si,])
      }
      #### control of the bad variance components = 0
      #### if during the first iteration there was an updated value < -.02
      #### and was a variance (not covariance) component, we set such var.comp
      #### for that trait to zero with their respective covariance (i.e. dominance) 
      if(length(bad2)>0){ #once they are zero keep them zero
        #print("error")
        #which traits ...theta[bad2,] <- 1e-8# 0
        trt.to.exclude <- vcv[bad2,1] #only first column since the 2nd is the same(are variances)
        trackin <- unlist(lapply(varos2,function(x){dim(x)[1]})); yas<- list()
        for(i in 1:length(trackin)){if(i==1){yas[[i]] <- 1:trackin[i]}else{yas[[i]]<- (trackin[i-1]+1):sum(trackin[1:i])}}
        vc.to.exc <- which(unlist(lapply(yas,function(x,y){fff <- which(x%in%y); if(length(fff)>0){return(TRUE)}else{return(FALSE)}},y=bad2)))
        vc.to.exc <- vc.to.exc[which(vc.to.exc != (length(ZETA)+1))]
        # now we not only adjust that var.comp but all its covar.comp
        var.com[vc.to.exc] <-lapply(var.com[vc.to.exc],function(x,y){
          matorral <- unique(t(combn(c(1:dim(x)[1],1:dim(x)[1]),2)))
          usethem <- matorral[which((matorral[,1] %in% y) | (matorral[,2] %in% y)),]
          x[usethem] <- 0#1e-7
          return(x)},y=trt.to.exclude)
        varos <- lapply(var.com, function(x){aa <- upper.tri(x); diag(aa) <- TRUE
        babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
        return(x[babas])})
        var.comM <- as.matrix(unlist(varos)) # current values of vc
      }
      ### end of the control for bad variance components
      
      var.comp.ret <- lapply(var.com, function(x,y,z){(x*y)/z},y=or.var,z=sc.var)
      #print(var.comp.ret)
      
      if(draw){
        asDF <- lapply(var.com,function(x){
          if(dim(as.matrix(x))[1]>1){base <- upper.tri(x)
        diag(base) <- TRUE
        base2 <- which(base,arr.ind = TRUE)
        base2 <- base2[ order(base2[,1], base2[,2]), ]
        return(as.matrix(x[base2]))
        }else{return(as.matrix(x))}})
        names(asDF) <- c(names(ZETA),"Residual")
        #print(var.coms)
        #print(asDF)
        var.coms <- mapply(cbind, var.coms, asDF, SIMPLIFY=FALSE)
        ## plot loglikelihood
        metric <- length(ZETA)+2
        is.even <- function(x) x %% 2 == 0 
        if(!is.even(metric)){metric <- metric+1}
        axo <- best.layout(metric)
        layout(matrix(1:metric,axo[1],axo[2]))
        plot(logLL,type="l", bty="n",col="cadetblue", lwd=2, main="LogLikelihood Multivariate\nAverage Information", las=2, xaxt="n",xlab="Iteration")
        axis(1,at=1:100,labels = 1:100)
        ## plot variance components
        if(nvarcom<2){nono=2}else{nono=nvarcom}
        palo <- brewer.pal(nono+1, "Accent")
        limos <- unlist(lapply(var.coms,function(x){c(max(x),min(x))})) # limits to draw
        ylimo <- c(min(limos), 1)
        
        for(g in 1:length(var.coms)){ # for each var.comp
          plot(var.coms[[1]][1,],type="l", bty="n",col="white", lwd=2, ylim=ylimo, ylab="Scaled Var.Comp values", xlab="Iteration", las=2,xaxt="n",main=paste("Scaled",varosss[g], "\nVar.Comp values"))
          axis(1,at=1:100,labels = 1:100)
          for(f in 1:(dim(var.coms[[g]])[1])){ # a.t1, at12, at2, b...
            lines(var.coms[[g]][f,],col=palo[g], lwd=2, lty=f) 
          }
          legend("topleft",legend=paste("t",traitm[,1],"t",traitm[,2],sep=""),lty=1:(dim(traitm)[1]),lwd=1.5, bty="n", cex=0.6, title="Var-Covar")
        }
        
      }
      ### ================================================= ### 
      ### ================================================= ###  
      
      ### variance components in the original scale (not standarized)
    }
    ## ================= END ALGORITHM ==================== ##
    ## ================= END ALGORITHM ==================== ##
    ## ================= END ALGORITHM ==================== ##
    ## ================= END ALGORITHM ==================== ##
    ## ================= END ALGORITHM ==================== ##
  }else{
    var.comp.ret <- forced
  }
  
  
  ### ================================================= ### 
  ### ================================================= ### 
  layout(matrix(1,1,1))
  names(var.comp.ret) <- varosss
  vmi.factor <- max(sc.var/or.var)
  # NEW V MATRIX WITH NORMAL VALUES
  # emmreml solve(kronecker(ZKZt, Vgt) + kronecker(diag(n), Vet) + tolparinv * diag(d * n))
  listGs <- list()
  for(i in 1:(nvarcom+1)){
    if(i<=nvarcom){ # ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
      zkz <- ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
      listGs[[i]] <- kronecker(as.matrix(zkz),var.comp.ret[[i]])
    }else{
      zkz <- diag(dimos[1]) #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
      listGs[[i]] <- kronecker(as.matrix(zkz),var.comp.ret[[i]])
    }
  }
  vm <- matrix(0,dimos[1]*ts,dimos[1]*ts)
  for(l in 1:length(listGs)){
    vm <- vm + listGs[[l]]
  }; vm[1:5,1:5];dim(vm)
  
  vmi <- try(solve(vm), silent = TRUE)
  if (class(vmi) == "try-error") {
    diag(vm) <- diag(vm) + (rep(1e-06, dim(vm)[2]))
    Vsing <- TRUE
    vmi <- solve(vm)
  }
  if(!silent){
    setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  }
  # END NEW V MATRIX WITH NORMAL VALUES
  
  
  ######## AIC BIC
  AIC = as.vector((-2 * logL ) + ( 2 * dim(X)[2]))
  BIC = as.vector((-2 * logL ) + ( log(length(Y)) * dim(X)[2]))
  ######## BETA.HAT, VAR.B.HAT
  xvx <- crossprod(X, vmi %*% X)
  xvxi <- solve(xvx) # variance of fixed effects
  pm <- vmi - vmi %*% X %*% solve(xvx, crossprod(X, vmi))
  beta <- xvxi %*% crossprod(X, vmi %*% Y.or3) # (XVX)-XV-Y .... Y.or3 %*% X %*% solve(crossprod(X))#
  rownames(beta) <- namesY
  varBhat <- solve(crossprod(X, vmi %*% X)) # XV-X'
  ######## E.HAT
  XB <- X%*%beta
  xb <- matrix(XB, nrow = nrow(Y.or), byrow = FALSE); #in a dataframe
  
  ## the fitted values will be done using the original X
  if(is.null(Xh)){ ##^^^
    XX <- matrix(rep(1,nrow(Yh)))  ##^^^
  }else{XX <- Xh} ##^^^
  XX <- do.call("adiag1", rep(list(XX), ts)) ##^^^
  xb.or <- XX%*%beta ##^^^
  xb.or <- matrix(xb.or, nrow = nrow(Yh), byrow = FALSE); ##^^^
  
  ehat <- matrix(t(Y.or2) - ((beta) %*% t(X.or)), ncol = 1, byrow = FALSE) #Y.or3 - XB # residuals = Y - XB
  popo <- dim(Y.or)[1]
  residu <- matrix(ehat, nrow = nrow(Y.or), byrow = TRUE); colnames(residu) <- namesY
  
  ######## U.HAT, PEV, etc.
  varvecG <- list()
  for (k in 1:nvarcom) { # ZKsZ'=G 
    K <-  ZETA[[k]]$K #tcrossprod(ZETA[[k]]$Z %*% ZETA[[k]]$K, (ZETA[[k]]$Z))
    varvecG[[k]] <- kronecker(as.matrix(K), (var.comp.ret[[k]]))
    #print(varvecG[[k]][1:5,1:5])
  }# str(varvecG)
  
  #print(dim(t(Y.or2)));print(dim(beta));print(dim(X.or)) ... 
  HobsInve <- vmi %*%  ehat # V-(y-Xb)' ... nxn %*% linearized(txn)
  #print(ya[1:5]);  print(HobsInve[1:5,])
  u.hat <- list()
  var.u.hat <- list()
  pev.u.hat <- list()
  
  dido <- dim(Yh) ## ^^^
  Zul <- list() ## ^^^
  
  for (k in 1:nvarcom) { # G ZV-(y-Xb)
    lev.re <- dim(ZETA[[k]]$Z)[2] # levels of the random effect
    Zforvec <- as(kronecker(t(as.matrix(ZETA[[k]]$Z)),diag(ts)), Class="sparseMatrix") # Z'
    #dim(varvecG[[k]])
    #dim((Zforvec))
    #dim(HobsInve)
    #dim(ZKforvec)
    ZKforvec <- varvecG[[k]] %*% Zforvec # KZ'
    #u.hats are returned mixed because the form of the V matrix
    #provi <- varvecG[[k]] %*% Zforvec %*% HobsInve # u.hat = GZ'V-(y-Xb)
    provi <- ZKforvec %*% HobsInve # u.hat = GZ'V-(y-Xb)
    var.u <- ZKforvec %*% pm %*% t(ZKforvec) # var.u.hat = ZGPZ'G ... sigma^4 ZKP ZK
    pev.u <- varvecG[[k]] - var.u #PEV.u.hat = G - ZGPGZ'
    #print(provi)
    u <- matrix(provi, nrow = lev.re, byrow = TRUE); colnames(u) <- namesY
    indnames <- colnames(ZETA[[k]]$Z)
    if(!is.null(indnames)){
      rownames(u) <- indnames
    }
    
    Zforvec.or <- as(kronecker((as.matrix(Zh[[k]]$Z)),diag(ts)), Class="sparseMatrix") ## ^^^
    zuu <-  (Zforvec.or %*% provi) ## ^^^
    Zul[[k]] <- matrix(zuu, nrow = lev.re, byrow = TRUE); colnames(Zul[[k]]) <- namesY ## ^^^
    
    u.hat[[k]] <- u
    var.u.hat[[k]] <- var.u
    pev.u.hat[[k]] <- pev.u
  }#str(u.hat)
  names(u.hat) <- varosss2
  names(var.u.hat) <- varosss2
  names(pev.u.hat) <- varosss2
  ##### COND. RESIDUALS AND FITTED
  
  ##### fitted Zu
  Zu <- do.call("+",Zul)
  fitted.y <- (xb.or + Zu)
  ## not really a Y.or is the one without missing data
  cond.ehat <- Y.or - fitted.y[ytouse,] # Y - (XB-Zu)
  
  dado <- lapply(ZETA, function(x){dim(x$Z)})
  if(AIsing){
    cat("\nInverse of the Average information matrix was singular")
  }
  if(Vsing){
    cat("\nInverse of V matrix was singular")
  }
  
  var.comp.ret <- lapply(var.comp.ret,function(x){round(x,7)})
  return(list(var.comp=var.comp.ret, V.inv=vmi, u.hat = u.hat , u.hat=u.hat,
              Var.u.hat = (var.u.hat), beta.hat = beta, Var.beta.hat = varBhat, 
              PEV.u.hat = pev.u.hat, residuals=residu, cond.residuals=cond.ehat,
              LL=logL, AIC=AIC, BIC=BIC, X=X, dimos=dado,
              fitted.y=fitted.y, fitted.u=Zu, ZETA=ZETA,
              method="MAI"))
}


