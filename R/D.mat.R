D.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
      n.core=1,shrink=FALSE,return.imputed=FALSE){
  
  #X <- apply(gg,2,function(x){y <- x; y[which(is.na(x))] <- mean(x, na.rm=TRUE); return(y)}); gg2[1:5,1:5]
  
  ty <- apply(X, 2, function(x){length(table(x))})
  vv <- which(ty == 3)
  X2 <- X[,vv]# only good markers with heterozygote plants
  # now transform 0 to 1's
  X3 <- apply(X2,2,function(x){y <- x; y[which(x == 1 | x==-1)] <- 0; y[which(x == 0)] <- 1; return(y)})
  X4 <- A.mat(X3, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
              n.core=n.core,shrink=shrink,return.imputed=return.imputed)
  return(X4)
}