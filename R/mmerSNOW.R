mmerSNOW <- function(y, X=NULL, Z=NULL, W=NULL, R=NULL, method="NR", REML=TRUE, iters=20, draw=FALSE, init=NULL, n.PC=0, P3D=TRUE, models="additive", ploidy=2, min.MAF=0.05, silent=FALSE, family=NULL, constraint=TRUE, sherman=FALSE, EIGEND=FALSE, Fishers=FALSE, gss=TRUE, forced=NULL, full.rank=TRUE, map=NULL, fdr.level=0.05, manh.col=NULL, gwas.plots=TRUE, tolpar=1e-6, tolparinv=1e-6){
  
  if(length(which(is.na(y))) == length(y)){
    stop("Y contains only missing data. Please double check your data.",call. = FALSE)
  }
  #if(is.null(Zi)){
  #  stop()
  #}
  #Z <- Zi
  ## convert to family function
  #if(!is.null(family)){
  # mox <- glm(y ~ 1, family = family)
  #mox$family[[1]]
  #y <- mox$family$linkfun(mox$y)
  #}
  # BLUPs for marker effects or BLUPs for breeding values are always specified in Z
  # extra function
  # needs: insert R matrices in EM and EMMA, AI??, correct var.beta.hat of fixed effects and PEV
  ### make sure X matrix from user is full rank when we remove the missing data
  
  if(!is.null(W)){
    if(dim(as.data.frame(y))[1] != dim(W)[1]){
      stop("The response 'Y' needs to have the same number of 
      observations than number of rows in the W matrix",call. = FALSE)
    }
    wnames <- colnames(W)
    if(is.null(wnames)){
      colnames(W) <- paste("M",1:(dim(W)[2]),sep="")
    }
    
  }
  ## make sure X is full rank
  if (!is.null(X)) {
    not.NA <- which(!is.na(y))
    Xtest <- as.matrix(X[not.NA, ])
    ################
    q <- qr(Xtest)
    chas <- q$pivot[seq(q$rank)]
    if(length(chas) < dim(X)[2]){
      if(!silent){
        cat("\nYour X matrix was not full rank, deleting columns to achieve full rank\n")
      }
      X <- as.matrix(X[,chas])
    }else{
      X <- as.matrix(X)
    }
  }
  ###
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  ### fix the input of the user
  if(is.list(Z)){
    if(is.list(Z[[1]])){ ### -- if is a 2 level list -- ##
      provided <- lapply(Z, names)
      for(s in 1:length(provided)){ #for each random effect =============================
        provided2 <- names(Z[[s]])
        if(length(provided2) ==1){ #----the 's' random effect has one matrix only----
          if(provided2 == "K"){ #user only provided K
            #zz <- diag(length(y))#model.matrix(~rownames(Z[[s]][[1]]))
            zz <- diag(length(y))
            #colnames(zz) <- rownames(Z[[s]][[1]])
            Z[[s]] <- list(Z=zz, K=Z[[s]][[1]])
          }
          if(provided2 == "Z"){ # user only provided Z
            #kk <- diag(dim(Z[[s]][[1]])[2])
            kk <- diag(dim(Z[[s]][[1]])[2])
            attributes(kk)$diagon <- TRUE
            #rownames(kk) <- colnames(Z[[s]][[1]]); colnames(kk) <- rownames(kk)
            Z[[s]] <- list(Z=Z[[s]][[1]],K=kk) 
          }
        }else{ #----the 's' random effect has two matrices----
          dido<-lapply(Z[[s]], dim) # dimensions of Z and K
          condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
          # condition, column size on Z matches with a square matrix K
          if(!condi){
            cat(paste("ERROR! In the",s,"th random effect you have provided or created an incidence \nmatrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2]," matrix \nas a variance-covariance matrix. Please double check your matrices.")
            stop()
          }
        }#---------------------------------------------------------------------------
      } #for each random effect end =================================================
    }else{ # if is a one-level list !!!!!!!!!!!!!
      if(length(Z) == 1){ ## -- if the user only provided one matrix -- ##
        provided <- names(Z)
        if(provided == "K"){
          #zz <- diag(length(y))
          zz <- diag(length(y))
          Z <- list(Z=zz, K=Z[[1]])
        }
        if(provided == "Z"){
          #kk <- diag(dim(Z[[1]])[2])
          kk <- diag(dim(Z[[1]])[2])
          attributes(kk)$diagon <- TRUE
          #rownames(kk) <- colnames(Z[[1]]); colnames(kk) <- rownames(kk)
          Z <- list(Z=Z[[1]],K=kk) 
        }
      }else{ # there's 2 matrices in Z
        dido<-lapply(Z, dim) # dimensions of Z and K
        condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
        # condition, column size on Z matches with a square matrix K
        if(!condi){
          cat(paste("ERROR! In the",s,"th random effect you have provided or created an incidence \nmatrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2]," matrix \nas a variance-covariance matrix. Please double check your matrices.")
          stop()
        }else{Z=list(Z=Z)}
      }
    }
  }else{
    if(is.null(Z)){ # the user is not using the random part
      cat("Error. No random effects specified in the model. \nPlease use 'lm' or provide a diagonal matrix in Z\ni.e. Zu = list(A=list(Z=diag(length(y))))\n")
      stop()
    }else{
      #stop;
      cat("\nThe parameter 'Z' needs to be provided in a 2-level list structure. \n\nPlease see help typing ?mmer and look at the 'Arguments' section\n")
      cat("\nIf no random effects provided, the model will be fitted using the 'lm' function\n\n")
    }
  }
  ###**********************************
  ##### FIX RANDOM EFFECTS FROM THE BEGGINING TO MAKE SURE LEVELS OF Z AND K ARE ORDERED
  Z <- lapply(Z, function(x){
    #########
    uuuz <- as.character(colnames(x$Z))#levels(as.factor(colnames(x$Z))) # order of Z
    uuuk <- as.character(colnames(x$K))#levels(as.factor(colnames(x$K))) # order of K
    #only if both Z and K have names, otherwise don't do anything
    if((length(uuuz) >0 ) & (length(uuuk) > 0)){ 
      inte <- intersect(uuuz,uuuk)
      if(length(inte)==length(uuuz)){ # the names were the same in Z and K
        x$K <- x$K[uuuz,uuuz]
      }else{ # no intersection between z and k names
        cat(paste("Column names of Z and K for random effect are not the same. Make sure they are in the correct order.\n"))
      }
    }
    return(x)
  })
  ###**********************************
  ## impute missing data in incidence matrices with median of the vector
  #Z <- lapply(Z, function(x){
  #  bbb <- which(is.na(x[[1]]))
  #  if(length(bbb)>0){
  #    im <- x[[1]] 
  #    im <- apply(im,2,function(y){qq <- which(is.na(y)); if(length(qq)>0){y[qq] <- median(y,na.rm=TRUE)}; return(y)})
  #    z=list(Z=im,K=x[[2]])
  #  }else{
  #    z <- x
  #  }
  #  return(z)
  #})
  ######## WRITE A POEM
  #if(!silent){
  #  poe(sample(1:9,1))
  #}
  
  ###
  ### if we want to get initial values based on lmer model we set lmerHELP=TRUE
  ###
  #if(lmerHELP){
  #  ## check that there's no matrices with class "hdm"
  #  were.hdm <- unlist(lapply(Z, function(x){if(!is.null(attributes(x$Z)$hdm)){y <-TRUE}else{y <- FALSE};return(y)}))
  #  
  #  if(length(which(were.hdm))==0){ # no hdm matrices then proceed
  #    
  #    ## only use it if the model has reps otherwise is a waste of time!!!
  #    not.square.track <- (unlist(lapply(Z,function(x){!is.square.matrix(x$Z)})))
  #    if(length(which(not.square.track)) > 0){ # check that at least any Z matrix is not square
  #      mmm <- length(Z)
  #      dddd <- vector(mode="list", length = mmm)
  #      #print(str(dddd))
  #      for(i in 1:mmm){ # for each random effect
  #        if(not.square.track[i]){ # if is not square proceed
  #          Z1 <- Z[[i]][[1]]
  #          #print(apply(Z1,1,function(x){which(x==1)}))
  #          dddd[[i]] <- colnames(Z1)[unlist(apply(Z1,1,function(x){which(x==1)}))]#
  #        }else{
  #          Z1 <- Z[[i]][[1]]
  #          #print(str(Z1))
  #          dddd[[i]] <- paste("LL",1:dim(Z1)[2],sep=".")
  #        }
  #      }
  #      dado <- as.data.frame(do.call("cbind",dddd))
  #      head(dado)
  #      # this returns which random effects can be evaluated by lmer, 1:YES, 0:NO
  #      doo <- apply(dado,2,function(x){s1<- length(x);s2<- length(unique(x));if(s1==s2){y <- FALSE}else{y<-TRUE}})
  #      goood <- which(doo)
  #      baaad <- which(!doo)
  #      doo2 <- names(dado)[which(doo)]
  #      a <- paste("y~",paste(paste("(1|",doo2,")",sep=""),collapse = "+"))
  #      lmermodel <- lmer(as.formula(a), data=dado)
  #      vc <- as.data.frame(VarCorr(lmermodel)); rownames(vc) <- vc$grp
  #      if(length(baaad)>0){ # if there was terms not possible to evaluate create a matrix for them
  #        bad.mat <- matrix(0,ncol=5,nrow=length(baaad)); rownames(bad.mat) <- names(baaad); colnames(bad.mat) <- colnames(vc)
  #        vc2 <- rbind(vc,bad.mat)
  #      }else{ # all was evaluated corectly
  #        vc2 <- vc
  #      }
  #      ## are all K matrices diagonals? in that case you could force values instead of initiating
  #      #were.squares <- unlist(lapply(Z,function(x){if(class(x$K)[1] == "ddiMatrix"){y <- TRUE}else{y<-FALSE};return(y)}))
  #      were.squares <- unlist(lapply(Z, function(x){if(!is.null(attributes(x$K)$diagon)){y <-TRUE}else{y <- FALSE};return(y)}))
  #      
  #      conditionK <- length(which(were.squares)) == length(Z)
  #      if(conditionK & (length(baaad) == 0)){ # if K's are square matrices and lmer evaluated everything force
  #        forced <- vc2[c(colnames(dado),"Residual"),4]
  #      }else{ # if Ks were not squared or lmer couldn't evaluate all just initiate
  #        init <- vc2[c(colnames(dado),"Residual"),4]
  #        were.zeros <- which(init == 0)
  #        if(length(were.zeros)>0){ # do not provide a zero or EM algorithm can fail
  #          init[were.zeros] <- .001
  #        }
  #      }
  #      
  #    }
  #    
  #  }
  #}
  ########### END OF LMER HELP
  #if(length(which(were.hdm))>0){ # rebaptize hdm matrices aas normal matrices
  #Z <- lapply(Z, function(x){if(class(x$Z)=="hdm"){class(x$Z) <- "matrix"};return(x)})
  #}
  ######## LET USER KNOW IF HE USED THE RIGHT PLOIDY LEVEL
  #if(!is.null(W)){
  #  ploidy.detected <- max(unlist(apply(W,2,function(x){abs(min(x)) + abs(max(x))})))
  #  if(ploidy != ploidy.detected){
  #    cat(paste("Please check the ploidy level selected, data indicates you may have\na ploidy level of",ploidy.detected,"and you specified a ploidy level:",ploidy,"\n\n"))
  #  }
  #}
  # X is fixed effects due to environmental factors
  # Z is random effects due to marker effects or genotype effects
  # W is an additional fixed effect due to markers if we have both experimental design 
  #    and markers and want to be estimated separately
  ## ------------------------------
  ## ------------------------------
  # if GWAS
  if(is.list(Z)){
    if(!is.null(Z) & !is.null(W)){
      #### impute response for GWAS in case will be used for bagging
      y[which(is.na(y))] <- mean(y, na.rm=TRUE)
      ####
      misso <- which(is.na(y))
      if(length(misso)>0){y[misso] <- mean(y,na.rm=TRUE)}
      #Wcnames <- colnames(W); Wrnames <- rownames(W)
      if(is.null(colnames(W))){colnames(W) <- paste("M",1:dim(W)[2],sep="-")}
      W <- apply(W, 2, function(x){
        vv <- which(is.na(x)); 
        if(length(vv) > 0){
          mu <- mean(x, na.rm = TRUE); 
          x[vv] <- mu}#else{x<-x}
        return(x)
      }
      )
      #colnames(W) <- Wcnames;rownames(W) <- Wrnames
      if(!silent){
        if(is.null(forced)){
          cat("Estimating variance components\n")
        }else{
          cat("Forcing variance components\n") 
        } 
      }
      fixed <- which(unlist(lapply(Z, function(x){names(x)[1]})) == "X") # elements of Z that are FIXED
      random <- which(unlist(lapply(Z, function(x){names(x)[1]})) == "Z") # elements of Z that are RANDOM
      random2 <- which(names(Z) == "Z")
      ## if Q+K model wants to be implemented
      if (n.PC > 0) {
        KK <- A.mat(W, shrink = FALSE)
        eig.vec <- eigen(KK)$vectors
        # if no X matrix make an intercept
        if(is.null(X)){X <- as.matrix(rep(1,dim(KK)[1]))}
        # X2 <- make.full(cbind(X[not.miss, ], Z %*% eig.vec[ix.pheno, 1:n.PC])) endelman
        #ZZ.comp <- list()
        Zss <- lapply(Z, function(x){x[[1]]})
        Zssp <- as(do.call("cbind", Zss),Class="sparseMatrix")
        #for(i in 1:length(Z)){
        #  X <- cbind(X, Z[[i]][[1]] %*% (1 + as.matrix(eig.vec[,1:n.PC]))) #plus 1 to ensure X is positive definite
        #}
        X <- make.full(cbind(X, Zssp %*% (1 + as.matrix(eig.vec[,1:n.PC])) ))
        #X <- make.full(X)
      }
      ## EMMA, EM, AI
      if(length(random) > 1 & method == "EMMA" & !silent){
        cat("For multiple random effects methods; 'AI', 'NR' and 'EM' are usually faster \nthan the 'EMMA' algorithm. Feel free to compare methods.\n")
      }
      if(method == "EMMA"){
        if(length(random) > 0){# EMMA both in 2-level list provided
          #res <- EMMA(y=y, X=X, Z=Z[[random]][[1]], K=Z[[random]][[2]],  REML=REML,silent=silent, EIGEND=EIGEND) 
          res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
          res$method <- method
          #names(res$u.hat) <- names(Z)[1]
          ### if user provide names in the Z component
          if(!is.null(names(Z))){
            rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
          }
          ### end of adding names
        }else{
          #res <- EMMA(y=y, X=X, Z=Z[[1]], K=Z[[2]],REML=REML,silent=silent,EIGEND=EIGEND)
          Z <- list(Z)
          res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
          res$method <- method
          #names(res$u.hat) <- names(Z)[1]
          ### if user provide names in the Z component
          if(!is.null(names(Z))){
            rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
          }
          ### end of adding names
        }
      }else if(method == "EM"){
        res <- EM(y=y, X=X, ETA=Z, init=init, iters = iters, REML=REML, draw=draw, silent=silent, forced=forced)
        res$method <- method
      }else if(method == "AI"){
        #res <- AI(y=y, X=X, ZETA=Z, R=R, REML=REML, draw=draw, silent=silent, iters = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, EIGEND=EIGEND,Fishers=Fishers)
        if(length(Z) == 1){ # if only one variance component, make sure is not Z and K both diags
          dias <- unlist(lapply(Z[[1]], function(x){if(dim(x)[1] == dim(x)[2]){y <- is.diagonal.matrix(x)}else{y <- FALSE};return(y)}))
          if(length(which(dias)) == 2){ # if K and Z are diagonals do EMMA
            #res <- EMMA(y=y, X=X, Z=Z[[random]][[1]], K=Z[[random]][[2]],REML=REML,silent=silent,EIGEND=EIGEND)
            res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
            res$method <- "EMMA"
            #names(res$u.hat) <- names(Z)[1]
            ### if user provide names in the Z component
            if(!is.null(names(Z))){
              rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
            }
            ### end of adding names
          }else{
            res <- AI(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, iters = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, EIGEND=EIGEND,Fishers=Fishers, forced=forced)
            res$method <- method
          }
        }else{ # if multiple variance components
          res <- AI(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, iters = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, EIGEND=EIGEND,Fishers=Fishers,forced=forced)
          res$method <- method
        }
      }else if(method == "NR"){
        if(length(Z) == 1){ # if only one variance component, make sure is not Z and K both diags
          dias <- unlist(lapply(Z[[1]], function(x){if(dim(x)[1] == dim(x)[2]){y <- is.diagonal.matrix(x)}else{y <- FALSE};return(y)}))
          if(length(which(dias)) == 2){ # if K and Z are diagonals do EMMA
            #res <- EMMA(y=y, X=X, Z=Z[[random]][[1]], K=Z[[random]][[2]],REML=REML,silent=silent,EIGEND=EIGEND)
            res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
            res$method<- "EMMA"
            #names(res$u.hat) <- names(Z)[1]
            ### if user provide names in the Z component
            if(!is.null(names(Z))){
              rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
            }
            ### end of adding names
          }else{
            if(is.null(R)){
              res <- NR(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, maxcyc = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, Fishers=Fishers, forced=forced, EIGEND=EIGEND)
              res$method <- method
            }else{
              res <- NRR(y=y,X=X,Z=Z,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=iters,draw=draw,constraint = constraint)
              res$method<-"NRR"
            }
          }
        }else{ # if multiple variance components
          if(is.null(R)){
            res <- NR(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, maxcyc = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, Fishers=Fishers, forced=forced, EIGEND=EIGEND)
            res$method <- method
          }else{
            res <- NRR(y=y,X=X,Z=Z,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=iters,draw=draw,constraint = constraint)
            res$method<-"NRR"
          }
        }
      }else{
        stop("Unrecognized method. Please select one of the methods; 'NR', 'AI', 'EMMA' or 'EM'. ",call. = FALSE)
      }
      #estimate variance components using G
      ##########################################
      #NANA <- which(!is.na(y))
      #y <- y[NANA]
      #if(length(NANA) > 0){y[NANA] <- mean(y, na.rm=TRUE)}
      #if(n.PC > 0){X2 <- X}else{ X2 <- make.full(res$X)}
      X2 <- res$X
      min.MAF = min.MAF; n <- length(y)
      # score calculation uses the H- matrix from the mixedmodel based on the A or G matrix
      # and H = Z K Z' + lambda*I. 
      # then it calculates beta = [XH-X']- XH-y, the residuals: e = y - XB
      # then the V(e) =  eH-e/(n-p) = SSe/(n-p), 
      # Var-Cov(Beta) =  Var(e) * [XH-X']-
      # then the F statistic as Beta^2/Var(Beta),,, x = (n-p) / (n-p + 1 * F)
      # and finally the -log10 of beta distribution [-log10(pbeta(x, v2/2, v1/2))]
      if(!silent){
        cat("\nPerforming GWAS")
      }
      max.geno.freq= 1 - min.MAF
      #n <- dim(W)[1]
      #badMAF <- apply(W,2,function(x, max.geno.freq,n){geno.freq <- table(x)/n; if((max(geno.freq) <= max.geno.freq) &  (min(geno.freq) >= min.MAF)){y <- 1}else{y<-0};return(y)}, max.geno.freq=max.geno.freq, n=n)
      #goodM <- which(badMAF==1)
      #W <- W[,goodM]
      W.scores <- list(NA)
      # draw layout
      if(length(models) > 2){layout(matrix(1:4,2,2))}else{layout(matrix(1:length(models),1,length(models)))}
      # run the different models
      ### QQ function
      qq <- function(scores) {
        remove <- which(scores == 0)
        if (length(remove) > 0) {
          x <- sort(scores[-remove], decreasing = TRUE)
        }
        else {
          x <- sort(scores, decreasing = TRUE)
        }
        n <- length(x)
        unif.p <- -log10(ppoints(n))
        plot(unif.p, x, pch = 16, xlab=expression(paste("Expected ",-log[10],"(p.value)")),
             ylab = expression(paste("Observed ",-log[10],"(p.value)")), col=transp("cadetblue"), main="QQ-plot")
        lines(c(0, max(unif.p, na.rm=TRUE)), c(0, max(unif.p, na.rm=TRUE)), lty = 2, lwd=2, col="blue")
      }
      ### END QQ FUNCTION
      
      deviations <- apply(W,2,sd,na.rm=TRUE) # sd of markers
      dev.no0 <- which(deviations > 0) # markers that are not singlular
      W <- W[,dev.no0] # only good markers will be tested
      
      for(u in 1:length(models)){ #### GWAS for all models specified ####
        
        model <- models[u]
        if(!silent){cat(paste("\nRunning",model,"model"))}
        ZO <- diag(dim(W)[1])
        step2 <- score.calc(marks=colnames(W),y=y,Z=ZO,X=X2,K=res$K, ZZ= res$Z, M=W,Hinv=res$V.inv,ploidy=ploidy,model=model,min.MAF=min.MAF,max.geno.freq=max.geno.freq, silent=silent, P3D=P3D, method=method)
        W.scores[[u]] <- as.matrix(step2$score)
        rownames(W.scores[[u]]) <- colnames(W)
        
        ##################### -----------------------------------------
        #### PLOTS
        ##################### -----------------------------------------
        if(!is.null(map) & gwas.plots){ 
          ########### MAP PRESENT  and user wants plots ##################
          dd <- W.scores[[u]]#matrix(step2$score)
          ffr <- fdr(dd, fdr.level=fdr.level)$fdr.10
          #rownames(dd) <- colnames(W)
          ## make sure map doesn't have duplicated markers
          non.dup <- which(!duplicated(map$Locus))
          map2 <- map[non.dup,]
          rownames(map2) <- map2$Locus
          ##get marker in common between GWAS and map 
          intro <- intersect(rownames(map2),rownames(dd))
          choco <- which(colnames(map2) == "Chrom")
          if(length(intro) > 0 & length(choco) > 0){ ####$$$$$ MARKERS IN COMMON  $$$$$$$$#######
            ## map adjusted and log p.values adjusted
            map3 <- map2[intro,]
            dd2 <- as.matrix(dd[intro,])
            map3$p.val <- dd[intro,]
            ## make plot
            if(is.null(manh.col)){
              col.scheme <- rep((transp(c("cadetblue","red"))),30)#heat.colors(12)#brewer.pal(12,"Accent")#
            }else{
              col.scheme <- rep(manh.col,30)#heat.colors(12)#brewer.pal(12,"Accent")#
            }
            layout(matrix(c(1,2,2),1,3))
            
            if(gwas.plots){ # user gave map, wants map, BUT WANTS PLOT??
              qq(step2$score)
              
              yylim <- ceiling(max(dd2,na.rm=TRUE))
              
              plot(dd2, bty="n", col=col.scheme[factor(map3$Chrom, levels = unique(map3$Chrom, na.rm=TRUE))], xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), pch=20, cex=2.5, las=2, ylim = c(0,yylim))
              ## make axis
              init.mrks <- apply(data.frame(unique(map3$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=map3$Chrom)
              fin.mrks <- apply(data.frame(unique(map3$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=map3$Chrom)
              inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
              
              axis(side=1, at=inter.mrks, labels=paste("Chr",unique(map3$Chrom),sep=""), cex.axis=.5)
              abline(h=ffr, col="slateblue4", lty=3, lwd=2)
              legend("topright", legend=paste("FDR(",fdr.level,")=",round(ffr,2), sep=""), 
                     bty="n", lty=3, lwd=2, col="slateblue4", cex=0.8)
            }
            
          }else{####$$$$$ NO MARKERS IN COMMON EXIST $$$$$$$$#######
            cat("\nError found! There was no markers in common between the column names of the W matrix \nand the map you provided. Please make sure that your data frame has names \n'Chrom' and 'Locus' to match correctly your map and markers tested. Plotting all markers.\n")
            map3 <- NULL
            layout(matrix(1:2,1,2))
            ffr <- fdr(step2$score, fdr.level=fdr.level)$fdr.10
            qq(step2$score)
            yylim <- ceiling(max(step2$score,na.rm=TRUE))
            plot(step2$score, col=transp("cadetblue", 0.6), pch=20, xlab="Marker index", 
                 ylab=expression(paste(-log[10],"(p.value)")), main=paste(model,"model"), bty="n", cex=1.5, ylim=c(0,yylim))
            abline(h=ffr, col="slateblue4", lty=3, lwd=2)
          }
          
        }else if (is.null(map) & gwas.plots){ 
          ############ NO MAP PROVIDED  but user wants plots#############
          layout(matrix(c(1,2,2),1,3))
          ffr <- fdr(step2$score, fdr.level=fdr.level)$fdr.10
          #layout(matrix(1:2,1,2))
          qq(step2$score)
          map3 <-NULL
          yylim <- ceiling(max(step2$score,na.rm=TRUE))
          plot(step2$score, col=transp("cadetblue", 0.6), pch=20, xlab="Marker index", 
               ylab=expression(paste(-log[10],"(p.value)")), main=paste(model,"model"), bty="n", cex=1.5)
          abline(h=ffr, col="slateblue4", lty=3, lwd=2)
        }else if(!is.null(map) & !gwas.plots){
          ## there IS map but user don't want plots
          dd <- W.scores[[u]]#matrix(step2$score)
          ffr <- fdr(dd, fdr.level=fdr.level)$fdr.10
          #rownames(dd) <- colnames(W)
          ## make sure map doesn't have duplicated markers
          non.dup <- which(!duplicated(map$Locus))
          map2 <- map[non.dup,]
          rownames(map2) <- map2$Locus
          ##get marker in common between GWAS and map 
          intro <- intersect(rownames(map2),rownames(dd))
          choco <- which(colnames(map2) == "Chrom")
          if(length(intro) > 0 & length(choco) > 0){ ####$$$$$ MARKERS IN COMMON  $$$$$$$$#######
            ## map adjusted and log p.values adjusted
            map3 <- map2[intro,]
            dd2 <- as.matrix(dd[intro,])
            map3$p.val <- dd[intro,]
          }else{####$$$$$ NO MARKERS IN COMMON EXIST $$$$$$$$#######
            cat("\nError found! There was no markers in common between the column names of the W matrix \nand the map you provided. Please make sure that your data frame has names \n'Chrom' and 'Locus' to match correctly your map and markers tested. Plotting all markers.\n")
            map3 <- NULL
          }
        }else if(is.null(map) & !gwas.plots){
          ## there is NO map and user don want plots
          map3 <- NULL
        }
        ####
        
      }
      names(W.scores) <- models
      #cat("\nMarker(s):\n")
      #anounce <- which(badMAF==0)
      #if(length(anounce)>0){print(colnames(W)[anounce])}
      #cat(paste("have been discarded. MAF was too low according to the threshold specified:",min.MAF))
      res$W.scores <- W.scores
      res$W <- W
      if(!is.null(map)){
        res$map <- map3
      }
      
    } # end of GWAS 
  }
  ## -------------------------------
  ## -------------------------------
  # if GENOMIC PREDICTION
  if(is.list(Z)){
    if((!is.null(Z) & is.null(W)) ) {
      
      if(!silent){
        if(is.null(forced)){
          cat("Estimating variance components\n")
        }else{
          cat("Forcing variance components\n") 
        } 
      }
      fixed <- which(unlist(lapply(Z, function(x){names(x)[1]})) == "X") # elements of Z that are FIXED
      random <- which(unlist(lapply(Z, function(x){names(x)[1]})) == "Z") # elements of Z that are RANDOM
      random2 <- which(names(Z) == "Z")
      ## EMMA, EM, AI
      if(length(random) > 1 & method == "EMMA" & !silent){
        cat("For multiple random effects methods; 'AI', 'NR' and 'EM' are usually faster \nthan the 'EMMA' algorithm. Feel free to compare methods.\n")
      }
      if(method == "EMMA"){
        if(length(random) > 0){# EMMA buth in 2-level list provided
          #res <- EMMA(y=y, X=X, Z=Z[[random]][[1]], K=Z[[random]][[2]],REML=REML, silent=silent,EIGEND=EIGEND) 
          res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
          res$method <- method
          #names(res$u.hat) <- names(Z)[1]
          ### if user provide names in the Z component
          if(!is.null(names(Z))){
            rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
          }
          ### end of adding names
        }else{
          #res <- EMMA(y=y, X=X, Z=Z[[1]], K=Z[[2]], REML=REML,silent=silent,EIGEND=EIGEND)
          Z <- list(Z)
          res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
          res$method <- method
          #names(res$u.hat) <- names(Z)[1]
          ### if user provide names in the Z component
          if(!is.null(names(Z))){
            rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
          }
          ### end of adding names
        }
      }else if(method == "EM"){
        res <- EM(y=y, X=X, ETA=Z, init=init, iters = iters, REML=REML, draw=draw, silent=silent, forced=forced)
        res$method <- method
      }else if(method == "AI"){
        if(length(Z) == 1){ # if only one variance component, make sure is not Z and K both diags
          dias <- unlist(lapply(Z[[1]], function(x){if(dim(x)[1] == dim(x)[2]){y <- is.diagonal.matrix(x)}else{y <- FALSE};return(y)}))
          if(length(which(dias)) == 2){ # if K and Z are diagonals do EMMA
            #res <- EMMA(y=y, X=X, Z=Z[[random]][[1]], K=Z[[random]][[2]],REML=REML,silent=silent,EIGEND=EIGEND)
            res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
            res$method <- "EMMA"
            #names(res$u.hat) <- names(Z)[1]
            ### if user provide names in the Z component
            if(!is.null(names(Z))){
              rownames(res$var.comp) <- c(paste("Var(",names(Z),")",sep=""), "Var(Error)")
            }
            ### end of adding names
          }else{
            res <- AI(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, iters = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, EIGEND=EIGEND,Fishers=Fishers, gss=gss, forced=forced)
            res$method <- method
          }
        }else{ # if multiple variance components
          res <- AI(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, iters = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, EIGEND=EIGEND,Fishers=Fishers, gss=gss, forced=forced)
          res$method <- method
        }
      }else if(method == "NR"){
        if(length(Z) == 1){ # if only one variance component, make sure is not Z and K both diags
          dias <- unlist(lapply(Z[[1]], function(x){if(dim(x)[1] == dim(x)[2]){y <- is.diagonal.matrix(x)}else{y <- FALSE};return(y)}))
          if(length(which(dias)) == 2){ # if K and Z are diagonals do EMMA
            #res <- EMMA(y=y, X=X, Z=Z[[random]][[1]], K=Z[[random]][[2]],REML=REML,silent=silent,EIGEND=EIGEND)
            #print(str(Z))
            res <- EMMA(y=y, X=X, ZETA=Z,  REML=REML,silent=silent, EIGEND=EIGEND, che=FALSE) 
            res$method<- "EMMA"
            #names(res$u.hat) <- names(Z)[1]
          }else{
            if(is.null(R)){
              #print("yes")
              res <- NR(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, maxcyc = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, Fishers=Fishers, forced=forced, EIGEND=EIGEND)
              res$method <- method
            }else{
              res <- NRR(y=y,X=X,Z=Z,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=iters,draw=draw,constraint = constraint)
              res$method<-"NRR"
            }
          }
        }else{ # if multiple variance components
          if(is.null(R)){
            res <- NR(y=y, X=X, ZETA=Z, REML=REML, draw=draw, silent=silent, maxcyc = iters, constraint=constraint, init=init, sherman=sherman, che=FALSE, Fishers=Fishers, forced=forced, EIGEND=EIGEND)
            res$method <- method
          }else{
            #print(str(Z))
            res <- NRR(y=y,X=X,Z=Z,R=R,tolpar=tolpar,tolparinv=tolparinv,maxcyc=iters,draw=draw,constraint = constraint)
            res$method<-"NRR"
          }
        }
      }else{
        stop("Unrecognized method. Please select one of the methods; 'NR', 'AI', 'EMMA' or 'EM'. ",call. = FALSE)
      }
      #res$method <- method
      res$maxim <- REML
      res$W <- W
    }
  }
  ## -------------------------------
  ## -------------------------------
  # if NULL MODEL or ONLY FIXED EFFECTS
  if((is.null(X) & is.null(Z) & is.null(W)) ) {
    res <- lm(y~1)
    jkl <- c(23,18,9,20,20,5,14, NA,2,25,NA,7,9,15,22,1,14,14,25,NA,3,15,22,1,18,18,21,2,9,1,19)
    oh.yeah <- paste(letters[jkl],collapse = "")
  }
  if((!is.null(X) & is.null(Z) & is.null(W)) ) {
    res <- lm(y~ X-1)
  }
  class(res)<-c("mmer")
  #res$u.hat <- lapply(res$u.hat, function(x,y){rownames(x) <- colnames(y); return(x)}, y=Z)
  layout(matrix(1,1,1))
  #if(!is.null(beeping)){
  #  beep(sound = beeping, expr = NULL) 
  #}
  return(res)
}
