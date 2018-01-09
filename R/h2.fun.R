
h2.fun <- function(object, data, gTerm=NULL, eTerm=NULL) {
  
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
    nn
    
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
    for(gs in geTerm){
      counter <- counter + 1
      if(counter==1){
        pev.mat <- pev.mat + object$PEV.u.hat[[gs]][[1]]
      }; pev.mat[1:4,1:4]
      A<- object$ZETA[[gs]]$K
      vc <- as.numeric(object$var.comp[[gs]])
      vcs[[gs]] <- vc
      G <- G + A * vc
    }
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
    M<-id.mat-(ginv%*%pev.mat)
    #M <- make.full(M)
    eM<-eigen(M)# eigenvalues of M
    sm<-rep(NA,nv)
    sm<-ifelse((Im(eM$values) ==0 |Re(eM$values) ==0 )==T, 1, 0)#number of non-zero eigenvectors
    neM<-sm*Re(eM$values) #  eigen values get a zero, and non-zero are multiplied by 1
    seM<-sum(neM) # add all eigen values, full heritability
    h2<- seM/sum(sm)
    
    h2.cullis <- 1 - (mean(diag(pev.mat))/vc)
    h2.stdrd <- as.numeric(vc/(vc+(ve/nn)))
    h2s[[e]] <- data.frame(H2.stdrd=h2.stdrd,H2.cullis=h2.cullis,H2.oakey.eigen=h2)
    
  }
  h2d <- as.data.frame(do.call(rbind,h2s))
  h2d[,eTerm] <- names(h2s)
  h2d[,"trait"] <- colnames(object$Y)
  return(h2d)
  
}
