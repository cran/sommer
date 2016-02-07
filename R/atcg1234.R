atcg1234 <- function (data, ploidy = 2, format = "ATCG", maf = 0) {
  apply_pb <- function(X, MARGIN, FUN, ...) {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    wrapper <- function(...) {
      curVal <- get("counter", envir = env)
      assign("counter", curVal + 1, envir = env)
      setTxtProgressBar(get("pb", envir = env), curVal + 
                          1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  user.code <- apply(data[, c(1:(round(dim(data)[2]/20)))], 
                     2, function(x) {
                       q <- which(!is.na(x))[1]
                       ss1 <- substr(x[q], start = 1, stop = 1)
                       ss2 <- substr(x[q], start = 2, stop = 2)
                       vv1 <- which(c(ss1, ss2) == "")
                       if (length(vv1) > 0) {
                         y <- 1
                       }
                       else {
                         y <- 0
                       }
                       return(y)
                     })
  AA <- sum(user.code, na.rm = TRUE)/length(user.code)
  if (AA > 0.9) {
    rrn <- rownames(data)
    gbs.to.bisnp <- function(x) {
      y <- rep(NA, length(x))
      y[which(x == "A")] <- "AA"
      y[which(x == "T")] <- "TT"
      y[which(x == "C")] <- "CC"
      y[which(x == "G")] <- "GG"
      y[which(x == "R")] <- "AG"
      y[which(x == "Y")] <- "CT"
      y[which(x == "S")] <- "CG"
      y[which(x == "W")] <- "AT"
      y[which(x == "K")] <- "GT"
      y[which(x == "M")] <- "AC"
      y[which(x == "+")] <- "NN"
      y[which(x == "0")] <- "NN"
      y[which(x == "-")] <- NA
      y[which(x == "N")] <- NA
      return(y)
    }
    cat("Converting GBS single-letter to biallelic code\n")
    data <- apply_pb(data, 2, gbs.to.bisnp)
    rownames(data) <- rrn
    data <- as.data.frame(data)
  }
  s1 <- rownames(data)
  s2 <- colnames(data)
  data <- as.data.frame(t(data))
  rownames(data) <- s2
  colnames(data) <- s1
  bases <- c("A", "C", "G", "T")
  get.ref <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ATCG") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        stop("Error in genotype matrix: More than 2 alleles")
      }
      if (sum(ans) == 2) {
        ref.alt <- bases[which(ans == 1)]
      }
      if (sum(ans) == 1) {
        ref.alt <- c(bases[which(ans == 1)], NA)
      }
    }
    return(ref.alt)
  }
  markers <- as.matrix(data)
  cat("Obtaining reference alleles\n")
  tmp <- apply_pb(markers, 1, get.ref, format = format)
  Ref <- tmp[1, ]
  Alt <- tmp[2, ]
  cat("Converting to numeric format\n")
  M <- apply_pb(cbind(Ref, markers), 1, function(x) {
    y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
    ans <- as.integer(lapply(y, function(z) {
      ifelse(z[1] < 0, ploidy, ploidy - length(z))
    }))
    return(ans)
  })
  gid.geno <- s1
  rownames(M) <- gid.geno
  bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
  if (bad > 0) {
    stop("Invalid marker calls.")
  }
  cat("Calculating minor allele frequency (MAF)\n")
  MAF <- apply_pb(M, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  polymorphic <- which(MAF > maf)
  M <- M[, polymorphic]
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  missing <- which(is.na(M))
  if (length(missing) > 0) {
    cat("Imputing missing data with mode \n")
    M <- apply_pb(M, 2, impute.mode)
  }
  if (ploidy == 2) {
    M <- M - 1
  }
  return(M)
}