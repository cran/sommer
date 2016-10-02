MNR <- function(Y,X=NULL,ZETA=NULL,R=NULL,init=NULL,maxcyc=20,tol=1e-3,tolparinv=1e-6,draw=TRUE,silent=FALSE, constraint=TRUE, EIGEND=FALSE, forced=NULL){
  
  choco <- names(ZETA)
  if(tol==1988.0906){
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
  #####((((((((((((((((((((((((((((((((()))))))))))))))))))))))))))))))))
  #namesX <- names(X)
  if(EIGEND){
    if(length(ZETA)>1){
      stop("The eigen decomposition to accelarate inversion only works for models with one relationship matrix",call. = FALSE)
    }
    casco <- dim(ZETA[[1]]$Z)
    if(casco[1] != casco[2]){
      stop("The eigen decomposition only works for square models",call. = FALSE)
    }
    cat("EIGEND feature activated. Eigen decomposition of K will be performed\n")
  }
  namesY <- colnames(Y)
  if(is.null(namesY)){
    namesY <- paste("T",1:dim(as.matrix(Y))[2],sep="")
  }
  
  Y <-apply((as.matrix(Y)),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  ts <- dim(as.matrix(Y))[2]
  base.var <- var(as.matrix(Y),na.rm=TRUE)
  sc.var <- var(as.matrix(scale(Y)),na.rm=TRUE)
  
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
  
  Y.or <- Y #apply((as.matrix(as.matrix(Y))),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  Y <- scale(Y)
  sonso <- Y
  
  if(is.null(X)){
    X <- matrix(rep(1,dim(as.data.frame(Y))[1])) 
  }
  
  if(EIGEND){
    EIGENS <- lapply(ZETA, function(x){eigen(x[[2]])}) # eigen decomposition of K
    Us <- lapply(EIGENS, function(x){x$vectors}) # extract eigen vectors U
    Usp <- as(do.call("adiag1", Us),Class="sparseMatrix") # U'G as diagonal
    Ds <- lapply(EIGENS, function(x){diag(x$values)}) # extract eigen values D
    Dsp <- as(do.call("adiag1", Ds),Class="sparseMatrix") # U'G as diagonal
    ZETA <- lapply(as.list(1:length(ZETA)),function(x,zz,kk){list(Z=zz[[x]][[1]], K=kk[[x]])}, zz=ZETA, kk=Ds)
    Y <- as.matrix(t(Usp) %*% as.matrix(Y))
    X <- as.matrix(t(Usp) %*% X)
    X.or <- X
    Y.or <- Y
  }
  
  Y.or2 <- as.matrix(as.vector(Y.or)); dim(Y.or2)
  
  n <- dim(ZETA[[1]]$Z)[1]#no. of individuals
  nvarcom <- length(ZETA)
  dimos <- dim(as.matrix(Y))
  ts <- dimos[2] #no. of traits
  inds <- dimos[1] #no. of individuals
  
  
  X.or <-X
  ###############
  # LINEARIZE
  ###############
  # impute phenotypes
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
  if(EIGEND){
    ZETA <- lapply(ZETA, function(x){x$Z <- as(x$Z,Class="sparseMatrix");x$K <- as(x$K,Class="sparseMatrix"); return(x)})
  }
  ########### list of ZKs to be used for derivatives in posterior calc
  listKs <- list()
  for(i in 1:(nvarcom+1)){
    if(i<=nvarcom){ # ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
      listKs[[i]] <- ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
    }else{
      if(!is.null(R)){
        rrr <- do.call("kronecker",R)#ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
        if(dim(rrr)[1] != dimos[1]){
          stop("Please check the residual structure passed to the function. When doing the\nkronecker product we did not find dimensions of R match with dimensions of Y.",call. = FALSE)
        }
        listKs[[i]] <- rrr #do.call("kronecker",R)#ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
      }else{
        listKs[[i]] <- diag(dimos[1]) #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) # 
      }
    }
  }
  ####################
  var.comp.ret<- list()
  
  if(!silent){
    count <- 0
    tot <- maxcyc
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  #cycle=4
  if(is.null(forced)){
    
    for(cycle in 1:maxcyc){
      
      if(!silent){
        count <- count + 1
      }
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
      
      V <- try(solve(as(W, Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
      if(class(V) == "try-error"){
        V <- try(solve(as(W + tolparinv * diag(dim(W)[2]), Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
      }
      dim(V) #dimens
      #print(V[1:5,1:5])
      W<-NULL
      
      # P = V-  -  V-X [X' V- X]-1 XV- =  WQK
      VX <- V %*% X # V- X
      P <- V - VX %*% solve(t(X)%*%VX, t(VX)) #WQK
      V <- NULL
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
      
      #sigma <- lapply(sigma,function(x){bad <- which(diag(x)<0); if(length(bad)>0){diag(x)[bad] <- 0.01};return(x)})
      #sigma[[1]][1,1:2]<-0
      #sigma[[1]][1:2,1]<-0
      
      var.comp.ret[[cycle]] <- lapply(sigma, function(x,y,z){(x*y)/z},y=base.var,z=sc.var)
      
      if(cycle > 1 & (delta.llik) < tol) {#tol*10
        if(!silent){
          setTxtProgressBar(pb, ((tot-1)/tot))### keep filling the progress bar
        }
        break
      }
      #cat(paste("LL=",llik,"iter=",cycle))
      llstore[cycle] <- llik
      ############ draw
      if(draw){
        layout(matrix(1:2,2,1))
        plot(llstore,type="l",lwd=2,bty="n",col="cadetblue",las=2,xaxt="n",main="LogLikelihood",ylab="LL", xlab="Iteration")
        axis(1,at=1:maxcyc,labels = 1:maxcyc)
        legend("topleft",legend=round(llik,3),bty="n",cex=0.8)
        plot(sigma3[1,],type="l",lwd=2,bty="n",ylim=c(min(sigma3),max(sigma3)),col="white",las=2,xaxt="n",main="Var-Covar components",ylab="scaled var.comp",xlab="Iteration")
        for(u in 1:dim(sigma3)[1]){
          lines(sigma3[u,],col="red",lty=u)
        }
        axis(1,at=1:maxcyc,labels = 1:maxcyc)
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
    sigma.scaled <- sigma
    good <- which(llstore==max(llstore))
    theta <- var.comp.ret[[good]]
    
    ############################################
    ############################################
    ############################################
    ############################################
    ############################################
    ## if there was negative variance components
    if(constraint){
      a <- lapply(theta,function(x){bad <- which(diag(x)<0,arr.ind = TRUE);return(bad)})
      where <- which(unlist(lapply(a, function(x){length(x)})) == 1)#tells you where is the bad random effect
      if(length(where) >0){
        if(!silent){
          setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
        }
        cat("\nOne or more variance components are zero. Boundary constraint applied\n")
        
        if(!silent){
          count <- 0
          tot <- maxcyc
          pb <- txtProgressBar(style = 3)
          setTxtProgressBar(pb, 0)
        }
        for(cycle in 1:maxcyc){
          
          for(u in 1:length(where)){
            uu <- where[u] # random effect
            to.fix <- a[[uu]] # trait
            for(m in 1:length(to.fix)){
              uuu <- to.fix[m]
              sigma[[uu]][uuu,] <- 0
              sigma[[uu]][,uuu] <- 0
            }
          }
          
          if(!silent){
            count <- count + 1
          }
          varos <- lapply(sigma, function(x){
            if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
            babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
            return(x[babas])}else{return(as.matrix(x))}})
          sigma2 <- as.matrix(unlist(varos)) # current values of vc
          
          #print(sigma)
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
          
          V <- try(solve(as(W, Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
          if(class(V) == "try-error"){
            V <- try(solve(as(W + tolparinv * diag(dim(W)[2]), Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
          }
          dim(V) #dimens
          V[1:5,1:5]
          W<-NULL
          
          # P = V-  -  V-X [X' V- X]-1 XV- =  WQK
          VX <- V %*% X # V- X
          P <- V - VX %*% solve(t(X)%*%VX, t(VX)) #WQK
          V <- NULL
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
          
          #sigma <- lapply(sigma,function(x){bad <- which(diag(x)<0); if(length(bad)>0){diag(x)[bad] <- 0.01};return(x)})
          #sigma[[1]][1,1:2]<-0
          #sigma[[1]][1:2,1]<-0
          
          var.comp.ret[[cycle]] <- lapply(sigma, function(x,y,z){(x*y)/z},y=base.var,z=sc.var)
          
          if(cycle > 1 & (delta.llik) < tol) {#tol*10
            ## force zeros last time
            for(u in 1:length(where)){
              uu <- where[u] # random effect
              to.fix <- a[[uu]] # trait
              for(m in 1:length(to.fix)){
                uuu <- to.fix[m]
                sigma[[uu]][uuu,] <- 0
                sigma[[uu]][,uuu] <- 0
              }
            }
            var.comp.ret[[cycle]] <- lapply(sigma, function(x,y,z){(x*y)/z},y=base.var,z=sc.var)
            
            if(!silent){
              setTxtProgressBar(pb, ((tot-1)/tot))### keep filling the progress bar
            }
            break
          }
          #cat(paste("LL=",llik,"iter=",cycle))
          llstore[cycle] <- llik
          ############ draw
          if(draw){
            layout(matrix(1:2,2,1))
            plot(llstore,type="l",lwd=2,bty="n",col="cadetblue",las=2,xaxt="n",main="LogLikelihood",ylab="LL", xlab="Iteration")
            axis(1,at=1:maxcyc,labels = 1:maxcyc)
            legend("topleft",legend=round(llik,3),bty="n",cex=0.8)
            plot(sigma3[1,],type="l",lwd=2,bty="n",ylim=c(min(sigma3),max(sigma3)),col="white",las=2,xaxt="n",main="Var-Covar components",ylab="scaled var.comp",xlab="Iteration")
            for(u in 1:dim(sigma3)[1]){
              lines(sigma3[u,],col="red",lty=u)
            }
            axis(1,at=1:maxcyc,labels = 1:maxcyc)
          }
          ########## end of draw
          
          ## end of one cycle
          if(!silent){
            setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
          }
          
        }
        #### end of cycles
        theta <- var.comp.ret[[cycle]]
        sigma.scaled <- sigma
      }
    }
    ############################################
    ############################################
    ############################################
    ############################################
    ############################################
    if(is.null(names(ZETA))){
      varosss <- c(paste("u.",1:length(ZETA), sep=""),"Residual")
      varosss2 <- c(paste("u.",1:length(ZETA), sep=""))
    }else{
      varosss <- c(names(ZETA),"Residual")
      varosss2 <- c(names(ZETA))
    }
    ############
    names(theta) <- varosss
    sigma <- lapply(theta, function(x,y){colnames(x)<-y;rownames(x)<-y;return(x)},y=namesY)
    
  }else{ #if sigma is forced
    sigma <- forced
    llik <- 1
    if(is.null(names(ZETA))){
      varosss <- c(paste("u.",1:length(ZETA), sep=""),"Residual")
      varosss2 <- c(paste("u.",1:length(ZETA), sep=""))
    }else{
      varosss <- c(names(ZETA),"Residual")
      varosss2 <- c(names(ZETA))
    }
  }
  ###################
  ## END OF VAR.COMP ESTIMATION
  ###################

  ###################
  ### 
  ### ================================================= ### 
  ### ================================================= ### 
  layout(matrix(1,1,1))
  listGs <- list()
  for(i in 1:(nvarcom+1)){
    if(i<=nvarcom){ # ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
      #zkz <- ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
      listGs[[i]] <- kronecker(as.matrix(listKs[[i]]),sigma[[i]])
    }else{
      #zkz <- diag(dimos[1]) #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z) #
      listGs[[i]] <- kronecker(as.matrix(listKs[[i]]),sigma[[i]])
    }
  }
  
  W <- matrix(0,dimos[1]*dimos[2],dimos[1]*dimos[2])
  for(l in 1:length(listGs)){ W <- W + listGs[[l]]};(W[1:5,1:5]);dim(W)
  
  vmi <- try(solve(as(W, Class="sparseMatrix"),sparse=TRUE), silent = TRUE)
  W <- NULL
  if(!silent){
    setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  }
  #dim(V) #dimens
  #vmi[1:5,1:5]
  #if(!silent){
  #  setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  #}
  # END NEW V MATRIX WITH NORMAL VALUES
  
  
  ######## AIC BIC
  AIC = as.vector((-2 * llik ) + ( 2 * dim(X)[2]))
  BIC = as.vector((-2 * llik ) + ( log(length(Y)) * dim(X)[2]))
  ######## BETA.HAT, VAR.B.HAT
  
  xvx <- crossprod(X, vmi %*% X)
  
  
  xvxi <- try(solve(xvx), silent = TRUE) # variance of fixed effects
  
  if(class(xvxi) == "try-error"){
    xvxi <- try(solve((xvx + tolparinv * diag(dim(xvx)[2]))), silent = TRUE)
  }
  
  pm <- vmi - vmi %*% X %*% solve(xvx, crossprod(X, vmi))
  beta <- xvxi %*% crossprod((X), vmi %*% Y.or2) # (XVX)-XV-Y .... Y.or3 %*% X %*% solve(crossprod(X))#

  #beta <- matrix(beta,nrow=nrow(beta))
  #rownames(beta) <- namesY
  #varBhat <- xvxi # XV-X'
  ######## E.HAT
  XB <- X%*%beta
  xb <- matrix(XB, nrow = nrow(Y.or), byrow = FALSE); #in a dataframe
  ##
  beta <- matrix(as.matrix(beta), nrow = ncol(X.or), byrow = FALSE); #in a dataframe
  #print(dim(t(Y.or)))
  #print(dim((t(beta) %*% t(X.or))))
  ##
  ehat <- matrix(t(Y.or) - (t(beta) %*% t(X.or)), ncol = 1, byrow = FALSE) #Y.or3 - XB # residuals = Y - XB
  residu <- matrix(ehat, nrow = nrow(Y.or), byrow = TRUE); colnames(residu) <- namesY
  ######## U.HAT, PEV, etc.
  varvecG <- list()
  for (k in 1:nvarcom) { # ZKsZ'=G 
    K <-  ZETA[[k]]$K #tcrossprod(ZETA[[k]]$Z %*% ZETA[[k]]$K, (ZETA[[k]]$Z))
    varvecG[[k]] <- kronecker(as.matrix(K), (sigma[[k]]))
    #print(varvecG[[k]][1:5,1:5])
  }# str(varvecG)
  #print(dim(t(Y.or2)));print(dim(beta));print(dim(X.or)) ... 
  HobsInve <- vmi %*%  ehat # V-(y-Xb)' ... nxn %*% linearized(txn)
  #print(ya[1:5]);  print(HobsInve[1:5,])
  u.hat <- list()
  var.u.hat <- list()
  pev.u.hat <- list()
  for (k in 1:nvarcom) { # GZ'V-(y-Xb)
    lev.re <- dim(ZETA[[k]]$Z)[2] # levels of the random effect

    Zforvec <- as(kronecker(t(as.matrix(ZETA[[k]]$Z)),diag(ts)), Class="sparseMatrix") # Z'

    ZKforvec <- varvecG[[k]] %*% Zforvec # GZ'
    #u.hats are returned mixed because the form of the V matrix
    
    if(EIGEND){ #kronecker(as.matrix(solve(t(Us[[k]]))),(diag(ts)))
      Uxi <- solve(t(Us[[k]]))# solve(t(Us[[k]]))
      provi <- (ZKforvec %*% HobsInve) # u.hat = GZ'V-(y-Xb) 
      u.hat[[k]] <- (Uxi) %*% matrix(provi, nrow = lev.re, byrow = TRUE); colnames(u.hat[[k]]) <- namesY
    }else{
      provi <- ZKforvec %*% HobsInve # u.hat = GZ'V-(y-Xb) 
      u.hat[[k]] <- matrix(provi, nrow = lev.re, byrow = TRUE); colnames(u.hat[[k]]) <- namesY
    }
    
    #print("a") #ZETA[[i]]$Z %*% tcrossprod(ZETA[[i]]$K, ZETA[[i]]$Z)
    #var.u.hat[[k]] <- ZKforvec %*% pm %*% t(ZKforvec) # var.u.hat = ZGPZ'G ... sigma^4 ZKP ZK
    var.u.hat[[k]] <- ZKforvec %*% tcrossprod(pm, ZKforvec) # var.u.hat = ZGPZ'G ... sigma^4 ZKP ZK
    
    #print("a")
    pev.u.hat[[k]] <- varvecG[[k]] - var.u.hat[[k]] #PEV.u.hat = G - ZGPGZ'
    #print("a")
    #print(provi)
    
    indnames <- colnames(ZETA[[k]]$Z)
    if(!is.null(indnames)){
      rownames(u.hat[[k]]) <- indnames
    }
  }#str(u.hat)
  names(u.hat) <- varosss2
  names(var.u.hat) <- varosss2
  names(pev.u.hat) <- varosss2
  ##### COND. RESIDUALS AND FITTED
  
  Zu <- matrix(0,dimos[1],dimos[2])
  for(o in 1:nvarcom){
    Zu <- Zu + ZETA[[o]]$Z %*% u.hat[[o]]
  }
  
  fitted.y <- (xb + Zu)
  cond.ehat <- Y.or - fitted.y # Y - (XB-Zu)
  dado <- lapply(ZETA, function(x){dim(x$Z)})
  
  sigma <- lapply(sigma,function(x){round(x,7)})
  
  ### fishers
  #sigma.cov <- (A.svd * 2)
  #print(A.svd * attr(sonso, 'scaled:scale') + attr(sonso, 'scaled:center'))
  if(is.null(forced)){
    sigma.cov <- ((A.svd) * 2) 
    FI <- (A)/2
    #### convert FI using pos
    FI.c <- matrix(0,dim(FI)[1],dim(FI)[2])
    FI.c <- FI / tcrossprod((sigmaxxx-1)*pos+1)
    ####print(A.svd)
    #####names(sigma) <- Vcoef.names
    sigma.cova <- try(ginv(FI.c),silent=TRUE)
  }else{
    sigma.cova <- NULL
  }
  
  
  colnames(beta) <- namesY
  #rownames(beta) <- namesX
  return(list(var.comp=sigma, V.inv=vmi, u.hat = u.hat , Var.u.hat = (var.u.hat), 
              beta.hat = beta, Var.beta.hat = xvxi, fish.inv=sigma.cova,
              PEV.u.hat = pev.u.hat, residuals=residu, cond.residuals=cond.ehat,
              LL=llik, AIC=AIC, BIC=BIC, X=X, dimos=dado, sigma.scaled=sigmaxxx,
              fitted.y=fitted.y, fitted.u=Zu, ZETA=ZETA,
              method="MNR",choco=choco))
}