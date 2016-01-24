fdr <- function(p, fdr.level=0.05, type="-log10") {
  ## this function takes p-values, and provides the -log10(p-val) adjusted for FDR
  ## use 10^-p to go back to p-values if you initially have -log10(p-val)
  smooth.df = 3
  if (min(p) < 0 || max(p) > 1) {
    print("ERROR: p-values not in valid range.")
    return(0)
  }
  lambda = seq(0, 0.9, 0.05)
  m <- length(p)
  pi0 <- rep(0, length(lambda))
  ## for each value of lambda, i.e. lambda=0.05
  ## mean(p > 0.05) / (1 - 0.05)
  for (i in 1:length(lambda)) {
    pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
  }
  spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
  pi0 <- predict(spi0, x = max(lambda))$y
  pi0 <- min(pi0, 1)
  if (pi0 <= 0) {
    print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
    return(0)
  }
  u <- order(p) # p values in order
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  qvalue[u[m]] <- min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                        1)
  }
  #### get the fdr in -log10 scale
  q.ans <- qvalue#(p=10^-input[, 4])
  if(type == "-log10"){
  temp <- cbind(q.ans, -log(p, base=10))#input[, 4]) # bind q.ans and -log10(p-values)
  }else{temp <- cbind(q.ans, p)}#input[, 4])}
  
  temp <- temp[order(temp[, 1]), ] # order them from smallest to largest
  
  temp2 <- tapply(temp[, 2], temp[, 1], mean) # for the "n" unique values in group1(q.vals) get the mean of p.values
  qvals <- as.numeric(rownames(temp2)) # which corresponds to the p-value <= to FDR
  x <- which.min(abs(qvals - fdr.level))
  first <- max(1, x - 2)
  last <- min(x + 2, length(qvals))
  if ((last - first) < 4) {
    last <- first + 3
  }
  splin <- smooth.spline(x = qvals[first:last], y = temp2[first:last], 
                         df = 3)
  res <- predict(splin, x = fdr.level)$y
  
  #lines(x = c(0, x.max), y = rep(predict(splin, x = fdr.level)$y, 2), lty = 2)
  
  return(res)
}
