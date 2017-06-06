pin <- function (object, transform){
  
  if(object$method %in% c("EMMA","EM")){
    stop("The pin function only works for 'NR' and 'AI' methods.",call. = FALSE)
  }
  
  if(object$method %in% c("MNR","MAI","MEMMA")){
    pframe <- as.list(object$sigma.scaled)#as.list(summary(object)[[3]][,1])
  }else{
    pframe <- as.list(object$var.comp[,1])
  }
 
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
                 pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  tname <- if (length(transform) == 3) 
    transform[[2]]
  else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$fish.inv
  toext <- upper.tri(Vmat)
  diag(toext) <- TRUE
  Vmat <- Vmat[which(toext,arr.ind = TRUE)]
  se <- sqrt(abs(sum(Vmat * X[i] * X[j] * k)))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}
