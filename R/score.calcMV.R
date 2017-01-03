score.calcMV <- function(marks,Y,Z,X,K,ZZ,M,Hinv,ploidy,model,min.MAF,max.geno.freq,silent=FALSE,P3D=TRUE) {
  # needs to be modified to work if Hinv is not provided, not in a hurry, normally wouldn't be possible
  #
  namesY <- names(Y)
  if(is.null(namesY)){
    namesY <- paste("T",1:dim(as.matrix(Y))[2],sep="")
  }
  ts<-dim(as.matrix(Y))[2]
  Y <-apply((as.matrix(Y)),2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return((x))})
  Y2 <- as.matrix(as.vector((Y))); dim(Y2)
  
  m <- length(marks)
  n <- nrow(Z)
  #ploidy=2
  #scores <- numeric(m)*NA
  v1 <- ncol(as.matrix(M[,marks[1]]))
  S <- as.matrix(M[,1])
  XX <- cbind(X,Z%*%S)
  BB <- ncol(XX)
  where.to.look <- dim(X)[2]+dim(Z%*%S)[2]
  where.betas <- seq(where.to.look,where.to.look*ts,where.to.look)
  scores <- matrix(NA,nrow=BB*ts,ncol=m)
  XX <- NULL
  beta.out <- scores
  #general <- length(grep("general",model,fixed=T))>0 
  ####################
  ## initialize the progress bar
  if(!silent){
    count <- 0
    tot <- m
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  #####################
  for (i in 1:m) {
    ###################
    if(!silent){
      count <- count + 1
    }
    ###################
    if(ploidy == 2){
      S <- as.matrix(M[,marks[i]])
    }else{
      S <- design.score(M[,marks[i]],model,ploidy,min.MAF,max.geno.freq)
    }
    
    if (!is.null(S)) {
      v1 <- ncol(S)
      X2 <- cbind(X,Z%*%S)
      p <- ncol(X2)
      v2 <- n - p  
      X2 <- do.call("adiag1", rep(list(X2), ts)); dim(X2)
      if (is.null(Hinv) | P3D == FALSE) {	# ZZ and K are only for y with no missing data
        yytt <- scale(Y)
        ETA <- list(list(Z=ZZ,K=K))
        out <- try(MNR(Y=yytt,X=X2,ZETA=ETA)) #Z=list(list(Z=ZZ,K=K)))
        #if (class(out)!="try-error") { 
        Hinv <- out$V.inv
        #}
      } 
      #print(str(X2))
      #print(str(Hinv))
      W <- crossprod(X2, Hinv %*% X2) 
      Winv <- try(solve(W),silent=TRUE)
      if (class(Winv) != "try-error") {
        beta <- Winv %*% crossprod(X2, Hinv %*% Y2)
        resid <- Y2 - X2 %*% beta
        s2 <- as.double(crossprod(resid, Hinv %*% resid))/v2
        Q <- s2 * Winv#[(p-v1+1):p,(p-v1+1):p]
        Q <- diag(diag(Q))
        Tt <- try(solve(Q), silent= TRUE)
        
        if (class(Tt) != "try-error") {
          V <- beta#[(p+1-v1):p]
          Fstat <- crossprod(t(diag(V[,1])),Tt%*%(diag(V[,1])))/v1
          x <- v2/diag(v2+v1*Fstat)
          #print(dim(scores))
          #print(length(x))
          scores[1:length(x),i] <- -log10(pbeta(x, v2/2, v1/2)) 
          beta.out[1:length(x),i] <- V[,1]
          #if (!general) {beta.out[i] <- beta[p]}                    
        }else{
          Tt <- try(solve(Q + (diag(dim(Q)[1])*1e-6) ), silent= TRUE)
          V <- beta#[(p+1-v1):p]
          Fstat <- crossprod(t(diag(V[,1])),Tt%*%(diag(V[,1])))/v1
          x <- v2/diag(v2+v1*Fstat)
          #print(dim(scores))
          #print(length(x))
          scores[1:length(x),i] <- -log10(pbeta(x, v2/2, v1/2)) 
          beta.out[1:length(x),i] <- V[,1]
          #if (!general) {beta.out[i] <- beta[p]}   
        }
      }
    }
    ################################
    if(!silent){
      setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    ################################
  }
  colnames(scores) <- marks
  scores <- scores[where.betas,]
  rownames(scores) <- namesY
  return(list(score=scores,beta=beta.out))			
}
