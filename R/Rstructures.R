## then can be multiplied by the error
AR1.mat = function(ar, dimo) {
  M = diag(dimo)
  M = ar^abs(row(M)-col(M))
  return(M)
}

ARMA.mat = function(ar, lam, dimo) {
  ## for ar
  M = diag(dimo)
  M = abs(row(M)-col(M))
  M[lower.tri(M)] <- M[lower.tri(M)]-1
  M[upper.tri(M)] <- M[upper.tri(M)]-1
  MM <- ar^M
  ## for lam
  N <- matrix(lam,dimo,dimo)
  diag(N) <- 0
  ## final
  MN <- MM*N
  return(MN)
}


ARMA.mat.dar = function(ar, lam, dimo) {
  ## for ar
  M = diag(dimo)
  M = abs(row(M)-col(M))
  M[lower.tri(M)] <- M[lower.tri(M)]-1
  M[upper.tri(M)] <- M[upper.tri(M)]-1
  MM <- M*ar^(M-1)
  ## for lam
  N <- matrix(lam,dimo,dimo)
  diag(N) <- 0
  ## final
  MN <- MM*N
  return(MN)
}

ARMA.mat.dlam = function(ar, lam, dimo) {
  ## for ar
  M = diag(dimo)
  M = abs(row(M)-col(M))
  M[lower.tri(M)] <- M[lower.tri(M)]-1
  M[upper.tri(M)] <- M[upper.tri(M)]-1
  MM <- ar^M
  ## for lam
  N <- matrix(1,dimo,dimo)
  diag(N) <- 0
  ## final
  MN <- MM*N
  return(MN)
}

# diagonal is like identity but each variance or diagonal element is different

# then can be multiplied by the error
CS.mat = function(ar, dimo) {
  M = matrix(ar,dimo,dimo)
  diag(M) <- 1
  return(M)
}

ID.mat = function(ar, dimo) {
  M = diag(dimo)
  return(M)
}
