score.calc <- function(marks,y,Z,X,K,ZZ,M,Hinv,ploidy,model,min.MAF,max.geno.freq,silent=FALSE,P3D=TRUE, method="NR") {
  # needs to be modified to work if Hinv is not provided, not in a hurry, normally wouldn't be possible
  #
  m <- length(marks)
  n <- nrow(Z)
  scores <- numeric(m)*NA
  beta.out <- scores
  general <- length(grep("general",model,fixed=T))>0 
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
      if (is.null(Hinv) | P3D == FALSE) {	# ZZ and K are only for y with no missing data
        yytt <- y[which(!is.na(y))]
        ETA <- list(list(Z=as.matrix(ZZ),K=as.matrix(K)))
        #print(str(ETA))
        #print(str(yytt))
        if(method=="EMMA"){
          out <- try(EMMA(y=yytt,X=X2,ZETA=ETA,silent=TRUE, che=FALSE)) #Z=list(list(Z=ZZ,K=K)))
        }else if(method=="NR"){
          out <- try(NR(y=yytt,X=X2,ZETA=ETA,silent=TRUE,che=FALSE, draw=FALSE)) 
        }else if(method=="AI"){
          out <- try(AI(y=yytt,X=X2,ZETA=ETA,silent=TRUE,che=FALSE, draw=FALSE)) 
        }else if(method=="EM"){
          out <- try(EM(y=yytt,X=X2,ZETA=ETA,silent=TRUE, draw=FALSE)) 
        }else{
          stop("Method not valid", call.=FALSE)
        }
                #print("yes")
        #if (class(out)!="try-error") {
        #print(str(out))
          Hinv <- out$V.inv
        #}
      } 
      W <- crossprod(X2, Hinv %*% X2) 
      Winv <- try(solve(W),silent=TRUE)
      if (class(Winv) != "try-error") {
        beta <- Winv %*% crossprod(X2, Hinv %*% y)
        resid <- y - X2 %*% beta
        s2 <- as.double(crossprod(resid, Hinv %*% resid))/v2
        Q <- s2 * Winv[(p-v1+1):p,(p-v1+1):p]
        Tt <- try(solve(Q),silent=TRUE)
        if (class(Tt) != "try-error") {
          V <- beta[(p+1-v1):p]
          Fstat <- crossprod(V,Tt%*%V)/v1
          x <- v2/(v2+v1*Fstat)
          mama <- -log10(pbeta(x, v2/2, v1/2)) 
          if(is.infinite(mama)){
            scores[i] <- NA
          }else{
            scores[i] <- mama
          }
          if (!general) {beta.out[i] <- beta[p]}                    
        }else{
          
            scores[i] <- NA
          
          if (!general) {beta.out[i] <- 1e-6} 
        }
      }
    }
    ################################
    if(!silent){
    setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    ################################
  }
  return(list(score=scores,beta=beta.out))			
}
