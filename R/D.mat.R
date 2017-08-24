D.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
                  n.core=1,shrink=FALSE,return.imputed=FALSE, ploidy=2, return.Xd=FALSE, 
                  method=1){
  
  #X <- apply(gg,2,function(x){y <- x; y[which(is.na(x))] <- mean(x, na.rm=TRUE); return(y)}); gg2[1:5,1:5]
  
  ty <- apply(X, 2, function(x){length(table(x))})
  vv <- which(ty > 2)
  
  if(length(vv)==0){
    cat("No heterozygous markers detected in the data. You might be using inbred lines.\nIf so, divide the markers in heterotic groups and do the kronecker product \namong A.mat's of the 2 groups to obtain a dominance relationship matrix\n")
    stop()
  }
  
  X2 <- X#[,vv]# only good markers with heterozygote plants
  # now transform 0 to 1's
  if(ploidy == 2){
    # Aa = 1 and AA|aa = 0
    X3 <- 1 - abs(X2)#apply(X2,2,function(x){y <- x; y[which(x == 1 | x==-1)] <- 0; y[which(x == 0)] <- 1; return(y)})
    if(return.Xd){
      X6 <- X3
    }else{
      if(method==1){
        X6 <- A.mat(X3, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                    n.core=n.core,shrink=shrink,return.imputed=return.imputed)
      }else if (method==2){
        
        M <- scale(X3, center = TRUE, scale = FALSE)
        K <- tcrossprod(M)
        
        bAlleleFrequency <- colMeans(X2+1)/2 # using original -1,0,1 matrix
        varHW <- sum((2 * bAlleleFrequency * (1 - bAlleleFrequency))^2) 
        
        X6 <- K/varHW
      }else if(method==3){
        print("using")
        #X3 <- 1 - abs(X2)
        n <- dim(X2)[1]
        p <- colSums(X2+1)/(2*n) # from marker marix in 0,1,2 format
        q <- 1-p
        varHW <- sum(2*p*q * (1-(2*p*q)) )
        
        X3pq <- X3 - (2*p*q)
        X6 <- tcrossprod(X3pq)/varHW
      }

    }
  }else{
    X3 <- X2 - (ploidy/2)
    possible <- (-(ploidy/2):(ploidy/2))
    homo <- c(possible[1],possible[length(possible)])
    hete <- setdiff(possible, homo)
    X4 <- apply(X3,2,function(x){y <- x; y[which(x %in% hete)] <- 1; y[which(x %in% homo)] <- 0; return(y)})
    ty2 <- apply(X4, 2, function(x){length(table(x))})
    vv2 <- which(ty2 > 1)
    X5 <- X4[,vv2]# only good
    if(return.Xd){
      X6 <- X5
    }else{
      if(method==1){
        X6 <- A.mat(X5, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                    n.core=n.core,shrink=shrink,return.imputed=return.imputed)
      }else if (method==2){
        
        M <- scale(X5, center = TRUE, scale = FALSE)
        K <- tcrossprod(M)
        
        bAlleleFrequency <- colMeans(X5+1)/2 # using original -1,0,1 matrix
        varHW <- sum((2 * bAlleleFrequency * (1 - bAlleleFrequency))^2) 
        
        X6 <- K/varHW
      }else if(method==3){
        #X3 <- 1 - abs(X2)
        n <- dim(X5)[1]
        p <- colSums(X5+1)/(2*n) # from marker marix in 0,1,2 format
        q <- 1-p
        varHW <- sum(2*p*q * (1-(2*p*q)) )
        
        X5pq <- X5 - (2*p*q)
        X6 <- tcrossprod(X5pq)/varHW
      }
      
    }
    
  }
  return(X6)
}