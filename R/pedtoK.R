pedtoK <- function(x, type="asreml"){
  
  if(type=="asreml"){
  K <- matrix(NA,max(x$Row),max(x$Column))
  for(i in 1:nrow(x)){
    K[x$Row[i],x$Column[i]] <- x$Ainverse[i]
  }
  K[1:3,1:3]
  copying <- function(m) { # copy upper triangular in lower triangular
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  copying2 <- function(m) { # copy lower triangular in upper triangular
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  
  K <- copying2(K)
  K[which(is.na(K), arr.ind = TRUE)] <- 0
  
  rownames(K) <- colnames(K) <- attr(x, "rowNames")
  
  Ks <- as(K, Class = "sparseMatrix")
  Ksi <- solve(Ks)
  rownames(Ksi) <- colnames(Ksi) <- attr(x, "rowNames")
  image(Ks)
  return(list(K=Ksi, Kinv=Ks))
  }
}