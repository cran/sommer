
h2.fun <- function(object, data, gTerm=NULL, eTerm=NULL, md=NULL) {
  
  if(missing(object)){
    stop("Please provide a model object of type mmer.\n", call. = FALSE)
  }
  if(is.null(gTerm)){
    stop("Please specify the gTerm in your model.\n", call. = FALSE)
  }
  if(missing(data)){
    stop("Please provide the dataset used to fit the model.\n", call. = FALSE)
  }
  
  if(!is.null(eTerm)){
    elevels <- as.character(na.omit(unique(data[,eTerm])))
  }else{
    data[,"FIELD"] <- "FIELD1"
    eTerm <- "FIELD"
    elevels <- "FIELD1"
  }
  
  h2s <- list()
  for(e in elevels){ # e <- elevels[1]
    if(length(elevels)==1){
      geTerm <- paste(gTerm,sep=":")
    }else{
      geTerm <- paste(e,gTerm,sep=":")
    }
    geTerm
    ## now get the heritability for each location
    v <- which(names(object$PEV.u.hat) %in% geTerm)
    if(length(v)==0){
      stop("The gTerm provided was not found in the model. You may not be providing the entire name.\n")
    }
    
    ## average sample size in such environment
    expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
    expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
    subdata <- droplevels(data[which(data[,eTerm]==e),])
    
    check <- grep("\\(",gTerm)
    if(length(check) > 0){
      gTerm2 <- expi(gTerm)
      nn <- median(table(subdata[,gTerm2]))
    }else{
      nn <- median(table(subdata[,gTerm]))
    }
    #print(nn)
    
    if(length(elevels)==1){ # single field
      f2 <- grep("units",names(object$var.comp))
      ve <- object$var.comp[[f2]]
    }else{
      f1 <- grep(e,names(object$var.comp))
      f2 <- grep("units",names(object$var.comp))
      ve <- object$var.comp[[intersect(f1,f2)]]
    }
    ve
    ##
    G <- 0
    pev.mat <- 0
    vcs <- numeric()
    counter <- 0
    meandiag <- numeric() # store the diagonal of the covariance matrix
    co <- 0
    for(gs in geTerm){
      co = co+1
      counter <- counter + 1
      if(counter==1){
        pev.mat <- pev.mat + object$PEV.u.hat[[gs]][[1]]
      }; pev.mat[1:4,1:4]
      A<- object$ZETA[[gs]]$K
      if(is.null(md)){
        meandiag[co] <- mean(diag(A))
      }else{meandiag[co] <- md}
      vc <- as.numeric(object$var.comp[[gs]])
      vcs[[gs]] <- vc
      G <- G + A * vc
    }
    meandiag <- mean(meandiag) # make sure we take the mean if more than one genetic term was fitted although that is not allowed

    G[1:3,1:3]
    try(ginv <- solve(G), silent = TRUE)
    if(class(ginv)=="try-error"){
      cat("Adding a small amount to the diagonal of A to make it positive-definite.\n")
      ginv <- solve(G+diag(1e-3,nrow(G)))
    }
    ginv[1:4,1:4]
    nv <- nrow(G)
    vc <- sum(as.numeric(vcs))
    
    id.mat<-diag(nv) #create identity matrix
    esh2<-1-(sum(pev.mat*ginv)/nv)
    #library(Matrix)
    M<-id.mat-(((1/meandiag)*ginv)%*%pev.mat)
    #M <- make.full(M)
    eM<-eigen(M)# eigenvalues of M
    sm<-rep(NA,nv)
    sm<-ifelse((Im(eM$values) ==0 |Re(eM$values) ==0 )==T, 1, 0)#number of non-zero eigenvectors
    neM<-sm*Re(eM$values) #  eigen values get a zero, and non-zero are multiplied by 1
    seM<-sum(neM) # add all eigen values, full heritability
    h2<- seM/sum(sm)
    
    pevs <- mean(diag(pev.mat))
    ve <- as.vector(ve)
    
    h2.cullis <- 1 - (mean(diag(pev.mat))/((meandiag)*vc))
    h2.stdrd <- as.numeric(vc/(vc+(ve/nn)))
    
    h2s[[e]] <- data.frame(PEV=pevs, Vg=vc, Ve=ve,N.rep=nn, H2.stdrd=h2.stdrd,H2.cullis=h2.cullis,H2.oakey.eigen=h2)
    # print(ve)
  }
  h2d <- as.data.frame(do.call(rbind,h2s))
  h2d[,eTerm] <- names(h2s)
  h2d[,"Trait"] <- colnames(object$Y)
  #h2d <- h2d[,c("Trait","Env",setdiff(colnames(h2d),c("Trait","Env")))]
  return(h2d)
  # pevs <- diag(mix2$PEV.u.hat$id$Yield)
  # vg <- mix2$var.comp$id[1,1]
  # ve <- mix2$var.comp$units[1,1]
  # pevm <- mix2$PEV.u.hat$id$Yield
  # m <- nrow(A)
  # Gi <- solve(diag(1,m)*vg)
  # 1 - mean(pevs/(2*vg)); mean(1 - pevs/(2*vg))
  # vg/(vg+ve)
  # uu <- (0.5*Gi %*% pevm)/m
  # 1 - matrix.trace(as.matrix(uu)); 1 - sum(diag(uu))
  # 1 - mean(diag(0.5*Gi %*% pevm))
}
