hdm <- function (data){
  v1 <- data$female
  v2 <- data$male
  factoro <- FALSE
  if ((!is.numeric(v1) | !is.numeric(v2))) {
    factoro <- TRUE
    vnames <- levels(as.factor(c(as.character(v1), as.character(v2))))
    v3 <- as.numeric(as.factor(c(as.character(v1), as.character(v2))))
    v4 <- as.factor(c(as.character(v1), as.character(v2)))
    v1.1 <- v3[1:length(v1)]
    v2.1 <- v3[(length(v1) + 1):length(v3)]
    v1 <- v1.1
    v2 <- v2.1
  }
  ncol <- max(c(v1, v2), na.rm = TRUE)
  nrow <- length(v1)
  nana <- sort(unique(c(v1, v2)), decreasing = FALSE)
  if (length(v1) != length(v2)) {
    stop
    print("your maternal and paternal vectors need to have the same size")
  }
  Z <- matrix(0, nrow = nrow, ncol = ncol)
  colnames(Z) <- nana
  for (i in 1:length(v1)) {
    ww1 <- which(colnames(Z) %in% v1[i])
    Z[i, ww1] <- 1
  }
  for (j in 1:length(v2)) {
    ww2 <- which(colnames(Z) %in% v2[j])
    Z[j, ww2] <- 1
  }
  if (factoro) {
    colnames(Z) <- vnames
  }
  return(Z)
}