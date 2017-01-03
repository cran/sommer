MMERM <- function (Y, X=NULL, Z=NULL, W=NULL, method="NR", tolpar = 1e-06, tolparinv = 1e-06, draw=TRUE, REML=TRUE, silent=FALSE, iters=15, init=NULL, che=TRUE, EIGEND=FALSE, forced=NULL, P3D=TRUE, models="additive", ploidy=2, min.MAF=0.05, gwas.plots=TRUE, map=NULL,manh.col=NULL,fdr.level=0.05, constraint=TRUE, IMP=TRUE, n.PC=0) {
  
#   print(str(Z))
#   print(str(X))
#   print(str(Y))
  ## for some reason I don't understand quite well the multivariate models can only estimate
  ## variance components in a steady way with V=kronecker(sigma,ZKZ) and dV/ds=kronecker(traits,ZKZ)
  ## but for getting the real V matrix is V=kronecker(ZKZ,sigma) which is not the same, why???? no idea yet.
  ## MAI2 and MNR are proof of that, took me a while to make them work.
  qq <- function(scores) {
    remove <- which(scores == 0)
    if (length(remove) > 0) {
      x <- sort(scores[-remove], decreasing = TRUE)
    }else {
      x <- sort(scores, decreasing = TRUE)
    }
    n <- length(x)
    unif.p <- -log10(ppoints(n))
    plot(unif.p, x, pch = 16, xlab=expression(paste("Expected ",-log[10],"(p.value)")),
         ylab = expression(paste("Observed ",-log[10],"(p.value)")), col=transp("cadetblue"), main="QQ-plot")
    lines(c(0, max(unif.p, na.rm=TRUE)), c(0, max(unif.p, na.rm=TRUE)), lty = 2, lwd=2, col="blue")
  }
  ########################
  if(is.list(Z)){
    if(is.list(Z[[1]])){ ### -- if is a 2 level list -- ##
      provided <- lapply(Z, names)
      for(s in 1:length(provided)){ #for each random effect =============================
        provided2 <- names(Z[[s]])
        if(length(provided2) ==1){ #----the 's' random effect has one matrix only----
          if(provided2 == "K"){ #user only provided K
            #zz <- diag(length(y))#model.matrix(~rownames(Z[[s]][[1]]))
            zz <- diag(dim(Y)[1])
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
          zz <- diag(dim(Y)[1])
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
  ########################  
  ########################
  ########################
  ### IF GENOMIC SELECTION
  ########################  
  ########################
  if(is.null(W)){
    
    if(method=="EMMA"){ ### EMMA
      if(!silent){
        if(length(Z)>1){
          cat("\nFor multivariate models with multiple random effect\nplease use method='NR' or 'AI' and argument MVM=TRUE.") 
        }
      }
      res <- MEMMA(Y, X=X, ZETA=Z, tolpar = tolpar, tolparinv = tolparinv, che=che)
    }else if(method=="AI"){ ### AI
      res <- MAI2(Y,X=X,ZETA=Z,init=init,maxcyc=iters,tol=tolpar*10000,draw=draw,silent=silent,constraint = constraint,EIGEND = EIGEND,forced=forced,IMP=IMP)
    }else if(method=="NR"){ ### NR tol=.001 = 1e-6*1000,,,,,,, tol=.01 = 1e-6*10000 
      res <- MNR(Y,X=X,ZETA=Z,init=init,maxcyc=iters,tol=tolpar*10000,draw=draw,silent=silent,constraint = constraint,EIGEND = EIGEND,forced=forced,IMP=IMP)
    }else if(method=="EM"){
      stop("Method 'EM' is not available currently for multivariate models.",call. = FALSE)
    }else{
      stop("Unrecognized method. Please select one of the methods; 'NR', 'AI', 'EMMA' or 'EM'. ",call. = FALSE)
    }
    ######################## 
  }else{
    ########################  
    ########################
    ### IF GWAS
    ########################  
    ########################
    
    ####-----------------------
    ####-----------------------
    ## if Q+K model wants to be implemented
    if (n.PC > 0) {
      
      make.full <- function(X) {
        svd.X <- svd(X)
        r <- max(which(svd.X$d > 1e-08))
        return(as.matrix(svd.X$u[, 1:r]))
      }
      
      KK <- A.mat(W, shrink = FALSE)
      eig.vec <- eigen(KK)$vectors
      # if no X matrix make an intercept
      if(is.null(X)){X <- as.matrix(rep(1,dim(KK)[1]))}
      # extract incidence matrices from random effects
      Zss <- lapply(Z, function(x){x[[1]]})
      Zssp <- as(do.call("cbind", Zss),Class="sparseMatrix")
      # bind X matrix provided by user, Z's of random effects, and eigen vectors
      X <- make.full(cbind(X, Zssp %*% (1 + as.matrix(eig.vec[,1:n.PC])) ))

    }
    ####-----------------------
    ####-----------------------
    if(method=="EMMA"){
      if(!silent){
        if(length(Z)>1){
          cat("\nFor multivariate models with multiple random effect\nplease use method='NR' or 'AI' and argument MVM=TRUE.") 
        }
      }
      res <- MEMMA(scale(Y), X=X, ZETA=Z, tolpar = tolpar, tolparinv = tolparinv, che=che)
    }else if(method=="AI"){
      res <- MAI2(scale(Y),X=X,ZETA=Z,init=init,maxcyc=iters,tol=tolpar*10000,draw=draw,silent=silent,constraint = constraint,EIGEND = EIGEND,forced=forced,IMP=IMP)
    }else if(method=="NR"){
      res <- MNR(scale(Y),X=X,ZETA=Z,init=init,maxcyc=iters,tol=tolpar*10000,draw=draw,silent=silent, constraint = constraint,EIGEND = EIGEND,forced=forced,IMP=IMP)
    }else if(method=="EM"){
      stop("Method 'EM' is not available currently for multivariate models.",call. = FALSE)
    }else{
      stop("Unrecognized method. Please select one of the methods; 'NR', 'AI', 'EMMA' or 'EM'. ",call. = FALSE)
    }
    ##################
    ### now get scores
    ##################
    deviations <- apply(W,2,sd,na.rm=TRUE) # sd of markers
    dev.no0 <- which(deviations > 0) # markers that are not singlular
    W <- W[,dev.no0] # only good markers will be tested
    
    rw <- rownames(W)
    cw <- colnames(W)
    W <- apply(W, 2, function(x){
      vv <- which(is.na(x)); 
      if(length(vv) > 0){
        mu <- mean(x, na.rm = TRUE); 
        x[vv] <- mu}#else{x<-x}
      return(x)
    }
    )
    rownames(W) <- rw
    colnames(W) <- cw
    
    if(is.null(colnames(W))){colnames(W) <- paste("M",1:dim(W)[2],sep="")}
    marks <-colnames(W)
    Y <- scale(Y)
    dimos <- dim(as.matrix(Y))
    Z <- diag(dimos[1])
    if(is.null(X)){X <- as.matrix(rep(1,dimos[1]))}
    K <- res$ZETA[[1]]$K
    ZZ <- res$ZETA[[1]]$Z
    M <- W
    Hinv <- res$V.inv
    max.geno.freq= 1 - min.MAF
    if(!silent){
      cat("\nPerforming GWAS\n")
    }
    W.scores <- score.calcMV(marks=marks,Y=Y,Z=Z,X=X,K=K,ZZ=ZZ,M=M,Hinv=Hinv,min.MAF=min.MAF,model=models[1],max.geno.freq=max.geno.freq,silent=silent,P3D=TRUE,ploidy=ploidy)
    res$W <- W
    res$W.scores <- W.scores
    ######################## 
    ########################
    #### plot
    ########################
    ########################
    if(!is.null(map) & gwas.plots){ ########### MAP PRESENT ##################
      dd <- t(W.scores$score)#matrix(W.scores$score)
      ffr <- apply(dd,2,function(x,y){fdr(x,fdr.level = y)$fdr.10},y=fdr.level)
      #ffr <- fdr(dd, fdr.level=fdr.level)$fdr.10
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
        
        for(t in 1:dim(W.scores$score)[1]){
          if(gwas.plots){ # user gave map, wants map, BUT WANTS PLOT??
            qq(W.scores$score[t,])
            yylim <- ceiling(max(dd2[,t],na.rm=TRUE))
            plot(dd2[,t], bty="n", main=colnames(dd2)[t], col=col.scheme[factor(map3$Chrom, levels = unique(map3$Chrom, na.rm=TRUE))], xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), pch=20, cex=2.5, las=2, ylim = c(0,yylim))
            ## make axis
            init.mrks <- apply(data.frame(unique(map3$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=map3$Chrom)
            fin.mrks <- apply(data.frame(unique(map3$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=map3$Chrom)
            inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
            
            axis(side=1, at=inter.mrks, labels=paste("Chr",unique(map3$Chrom),sep=""), cex.axis=.5)
            abline(h=ffr[t], col="slateblue4", lty=3, lwd=2)
            legend("topright", legend=paste("FDR(",fdr.level,")=",round(ffr[t],2), sep=""), 
                   bty="n", lty=3, lwd=2, col="slateblue4", cex=0.8)
          }
        }
        
        res$map <- cbind(map3[,1:3], map3$p.val)# map3
        
      }else{####$$$$$ NO MARKERS IN COMMON EXIST $$$$$$$$#######
        cat("\nError found! There was no markers in common between the column names of the W matrix \nand the map you provided. Please make sure that your data frame has names \n'Chrom' and 'Locus' to match correctly your map and markers tested. Plotting all markers.\n")
        map3 <- NA
        layout(matrix(1:2,1,2))
        for(t in 1:dim(W.scores$score)[1]){
          qq(W.scores$score[t,])
          yylim <- ceiling(max(W.scores$score[t,],na.rm=TRUE))
          plot(W.scores$score[t,], col=transp("cadetblue", 0.6), pch=20, xlab="Marker index", 
               ylab=expression(paste(-log[10],"(p.value)")), main=rownames(W.scores$score)[t], bty="n", cex=1.5, ylim=c(0,yylim))
          abline(h=ffr[t], col="slateblue4", lty=3, lwd=2)
          legend("topright", legend=paste("FDR(",fdr.level,")=",round(ffr[t],2), sep=""), 
                 bty="n", lty=3, lwd=2, col="slateblue4", cex=0.8)
        }
        res$map <- map3
      }
      
    } else if (is.null(map) & gwas.plots){ ############ NO MAP PROVIDED #############
      layout(matrix(c(1,2,2),1,3))
      dd <- t(W.scores$score)#matrix(W.scores$score)
      # get FDR value
      ffr <- apply(dd,2,function(x,yyy){fdr(x,fdr.level = yyy)$fdr.10},yyy=fdr.level)
      
      
      #layout(matrix(1:2,1,2))
      for(t in 1:dim(W.scores$score)[1]){
        qq(W.scores$score[t,])
        yylim <- ceiling(max(W.scores$score[t,],na.rm=TRUE))
        plot(W.scores$score[t,], col=transp("cadetblue", 0.6), pch=20, xlab="Marker index", 
             ylab=expression(paste(-log[10],"(p.value)")), main=rownames(W.scores$score)[t], bty="n", cex=1.5, ylim=c(0,yylim))
        abline(h=ffr[t], col="slateblue4", lty=3, lwd=2)
      }
      map3 <-NA
      res$map <- map3
      
    }
    ########################
    ########################  
    #### end of plot
    ########################
    ########################
  }
  
  return(res)
}