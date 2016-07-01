MMERM <- function (Y, X=NULL, Z=NULL, method="EMMA", tolpar = 1e-06, tolparinv = 1e-06, che=TRUE, silent=TRUE) {
  
  #if(method=="EMMA"){
  if(!silent){
    if(length(Z)>1){
      cat("Currently multivariate models are only available for a single random effect\nMerging all random effects in a single effect.") 
    }
  }
    res <- EMMAM(Y, X=X, ZETA=Z, tolpar = tolpar, tolparinv = tolparinv, che=che)
  #}
  return(res)
}