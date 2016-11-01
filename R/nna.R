nna<- function (pheno, trait = "y", rown = "row", coln = "col", nrows = 1, 
                ncols = 2) 
{
  bases <- c(rown, coln, trait)
  into <- intersect(colnames(pheno), bases)
  if (length(into) < 3) {
    stop("Please provide the information 'rown', 'coln' and 'y'.", 
         call. = FALSE)
  }
  newd <- pheno[, c(rown, coln, trait)]
  rn <- rownames(newd)
  newd <- apply(newd, 2, as.numeric)
  rownames(newd) <- rn
  newd <- data.frame(newd)
  pheno[, c(rown, coln, trait)] <- newd
  pheno$nnx <- NA
  nnx <- apply(newd, 1, function(x) {
    s1 <- as.numeric(x[1]-nrows)
    s2 <- as.numeric(x[1] - 1)
    s3 <- as.numeric(x[1] + 1)
    s4 <- as.numeric(x[1] + nrows)
    ns <- c(s1:s2, s3:s4)
    
    s5 <- as.numeric(x[2] - ncols)
    s6 <- as.numeric(x[2] - 1)
    s7 <- as.numeric(x[2] + 1)
    s8 <- as.numeric(x[2] + ncols)
    we <- c(s5:s6, s7:s8)
    
    
    ns <- ns[which(ns > 0)]
    we <- we[which(we > 0)]
    v1 <- which((newd[, rown] %in% ns) & (newd[, coln] %in% x[2]))
    v2 <- which(newd[, coln] %in% we & newd[, rown] %in% x[1])
    if (length(v1) > 0 & length(v2) > 0) {
      yy <- x[3] - apply(newd[c(v1, v2), ], 2, mean, na.rm = TRUE)[3]
    }  else if (length(v1) > 0 & length(v2) == 0) {
      yy <- x[3] - apply(newd[c(v1), ], 2, mean, na.rm = TRUE)[3]
    }  else if (length(v1) > 0 & length(v2) == 0) {
      yy <- x[3] - apply(newd[c(v2), ], 2, mean, na.rm = TRUE)[3]
    } else{
      yy<-x[3]
    }
    return(yy)
  })
  pheno$nnx <- as.numeric(unlist(nnx))
  return(pheno)
}