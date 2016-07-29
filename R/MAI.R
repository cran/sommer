MAI <- function(Y, X=NULL, ZETA=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=20, init=NULL, che=TRUE, EIGEND=FALSE, forced=NULL){
  
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
    varosss <- c(paste("u.",1:length(ZETA), sep=""),"Error")
    varosss2 <- c(paste("u.",1:length(ZETA), sep=""))
  }else{
    varosss <- c(names(ZETA),"Error")
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
      
      if (((logL-logL2<0.001) | (wi==iters)) ) {
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
        names(asDF) <- c(names(ZETA),"Error")
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
    #var.u <- matrix(provi2, nrow = lev.re, byrow = FALSE); #colnames(var.u) <- namesY
    #pev.u <- matrix(provi3, nrow = lev.re, byrow = FALSE); #colnames(pev.u) <- namesY
    
    u.hat[[k]] <- u
    var.u.hat[[k]] <- var.u
    pev.u.hat[[k]] <- pev.u
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
  cond.ehat <- Y.or2 - fitted.y # Y - (XB-Zu)
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


