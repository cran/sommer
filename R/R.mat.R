R.mat <- function(dims, rhor=.95, rhoc=.95, method="SIM"){
  car2 = function(dim, rhor, rhoc) {
    M = diag(dim)
    rhor^(row(M) - 1) * rhoc^(col(M) - 1)
  }
  car3 = function(dim, rhor, rhoc) {
    M = diag(dim)
    rhor^(row(M) - 1) #* rhoc^(col(M) - 1)
  }
  ######### symmetric matrix
  if(method=="SIM"){
    nn <- car2(dims,rhor, rhoc)
    
    dis <- (dim(nn)[1])
    s1 <- nn[1,]
    for(i in 2:dis){
      a1 <- s1[length(s1)]
      s1 <-c(a1,s1[-length(s1)])
      nn[i,] <- s1
    }
    nn2 <- t(nn)
    nn[lower.tri(nn)] <- nn2[lower.tri(nn2)]
  }
  ### autoregressive one direction
  if(method=="AR1"){
    nn <- car3(dims,rhor, rhoc)
  }
  ### autoregressive both directions
  if(method=="AR2"){
    nn <- car2(dims,rhor, rhoc)
  }
  return(nn)
}