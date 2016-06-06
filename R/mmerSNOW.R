mmerSNOW <-function (y, X = NULL, Z = NULL, W = NULL, R = NULL, method = "AI", 
          REML = TRUE, iters = 30, draw = FALSE, init = NULL, n.PC = 0, 
          P3D = TRUE, models = "additive", ploidy = 2, min.MAF = 0.05, 
          silent = FALSE, family = NULL, constraint = TRUE, sherman = FALSE, 
          EIGEND = FALSE, Fishers = FALSE, gss = TRUE, forced = NULL, 
          full.rank = TRUE, map = NULL, fdr.level = 0.05, manh.col = NULL, 
          gwas.plots = TRUE, lmerHELP = FALSE) 
{
  if (!is.null(X)) {
    case <- dim(X)[2]
    case2 <- dim(X)[2]
    xor <- X
    p <- ncol(X)
    not.NA <- which(!is.na(y))
    X <- as.matrix(X[not.NA, ])
    XtX <- crossprod(X, X)
    rank.X <- qr(XtX)$rank
    if (rank.X < p & full.rank) {
      while (rank.X < p) {
        case <- case - 1
        XtX <- crossprod(X[, -c(case:case2)], X[, -c(case:case2)])
        p <- ncol(X[, -c(case:case2)])
        rank.X <- qr(XtX)$rank
      }
      X <- as.matrix(xor[, -c(case:case2)])
      cat("\nYour X matrix is not full rank, deleting columns until full rank is achieved\n")
    }
    else {
      X <- xor
    }
  }
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  if (is.list(Z)) {
    if (is.list(Z[[1]])) {
      provided <- lapply(Z, names)
      for (s in 1:length(provided)) {
        provided2 <- names(Z[[s]])
        if (length(provided2) == 1) {
          if (provided2 == "K") {
            zz <- diag(length(y))
            Z[[s]] <- list(Z = zz, K = Z[[s]][[1]])
          }
          if (provided2 == "Z") {
            kk <- diag(dim(Z[[s]][[1]])[2])
            attributes(kk)$diagon <- TRUE
            Z[[s]] <- list(Z = Z[[s]][[1]], K = kk)
          }
        }
        else {
          dido <- lapply(Z[[s]], dim)
          condi <- (dido$Z[2] == dido$K[1] & dido$Z[2] == 
                      dido$K[2])
          if (!condi) {
            cat(paste("ERROR! In the", s, "th random effect you have provided or created an incidence \nmatrix with dimensions:", 
                      dido$Z[1], "rows and", dido$Z[2], "columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions", 
                      dido$Z[2], "x", dido$Z[2]), ", but you provided a", 
                dido$K[1], "x", dido$K[2], " matrix \nas a variance-covariance matrix. Please double check your matrices.")
            stop()
          }
        }
      }
    }
    else {
      if (length(Z) == 1) {
        provided <- names(Z)
        if (provided == "K") {
          zz <- diag(length(y))
          Z <- list(Z = zz, K = Z[[1]])
        }
        if (provided == "Z") {
          kk <- diag(dim(Z[[1]])[2])
          attributes(kk)$diagon <- TRUE
          Z <- list(Z = Z[[1]], K = kk)
        }
      }
      else {
        dido <- lapply(Z, dim)
        condi <- (dido$Z[2] == dido$K[1] & dido$Z[2] == 
                    dido$K[2])
        if (!condi) {
          cat(paste("ERROR! In the", s, "th random effect you have provided or created an incidence \nmatrix with dimensions:", 
                    dido$Z[1], "rows and", dido$Z[2], "columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions", 
                    dido$Z[2], "x", dido$Z[2]), ", but you provided a", 
              dido$K[1], "x", dido$K[2], " matrix \nas a variance-covariance matrix. Please double check your matrices.")
          stop()
        }
        else {
          Z = list(Z = Z)
        }
      }
    }
  }
  else {
    if (is.null(Z)) {
      cat("Error. No random effects specified in the model. \nPlease use 'lm' or provide a diagonal matrix in Z\ni.e. Zu = list(A=list(Z=diag(length(y))))\n")
      stop()
    }
    else {
      cat("\nThe parameter 'Z' needs to be provided in a 2-level list structure. \n\nPlease see help typing ?mmer and look at the 'Arguments' section\n")
      cat("\nIf no random effects provided, the model will be fitted using the 'lm' function\n\n")
    }
  }
  Z <- lapply(Z, function(x) {
    if (length(x) > 1) {
      if (length(x) > 1 | !is.diagonal.matrix(x[[2]])) {
        if (!is.null(colnames(x[[1]])) & !is.null(colnames(x[[2]]))) {
          if (length(which(colnames(x[[1]]) == colnames(x[[2]]))) != 
              dim(x[[1]])[2]) {
            y <- colnames(x[[1]])
            mo <- colnames(x[[2]])
            if (!is.null(mo)) {
              real <- apply(data.frame(mo), 1, function(b, 
                                                        z) {
                grep(b, z)[1]
              }, z = y)
              x[[1]] <- x[[1]][, real]
            }
            else {
              y2 <- strsplit(y, split = "")
              y3 <- y2[[1]]
              for (i in 2:length(y2)) {
                basd <- vector()
                for (j in 1:length(y3)) {
                  good <- which(y3[j] == y2[[i]])
                  if (length(good) > 0) {
                    basd[j] <- (y3[good])[1]
                    y2[[i]][good] <- NA
                  }
                }
                y3 <- basd
              }
              extraname <- paste(na.omit(y3), collapse = "")
              if (extraname != "") {
                real1 <- match(colnames(x[[2]]), gsub(as.character(extraname), 
                                                      "", as.character(colnames(x[[1]]))))
                if (length(which(is.na(real1))) == 0) {
                  x[[1]] <- x[[1]][, real1]
                }
                else {
                  cat("\nNames in Z and K matrices do not match. Might be that column names of your Z matrix \nare in different order than column names of K, or the column names of the Z matrix have\nan extra word. The analysis will be performed assuming column names of Z correspond \nto colum names of K and they are in the same order. Please take a look to make sure \nthe levels of K and Z are in the same order or do not have extra words.\n")
                }
              }
            }
          }
        }
      }
    }
    return(x)
  })
  Z <- lapply(Z, function(x) {
    bbb <- which(is.na(x[[1]]))
    if (length(bbb) > 0) {
      im <- x[[1]]
      im <- apply(im, 2, function(y) {
        qq <- which(is.na(y))
        if (length(qq) > 0) {
          y[qq] <- median(y, na.rm = TRUE)
        }
        return(y)
      })
      z = list(Z = im, K = x[[2]])
    }
    else {
      z <- x
    }
    return(z)
  })
  if (lmerHELP) {
    were.hdm <- unlist(lapply(Z, function(x) {
      if (!is.null(attributes(x$Z)$hdm)) {
        y <- TRUE
      }
      else {
        y <- FALSE
      }
      return(y)
    }))
    if (length(which(were.hdm)) == 0) {
      if (!is.square.matrix(Z[[1]]$Z)) {
        mmm <- length(Z)
        dddd <- vector(mode = "list", length = mmm)
        for (i in 1:mmm) {
          Z1 <- Z[[i]][[1]]
          dddd[[i]] <- colnames(Z1)[apply(Z1, 1, function(x) {
            which(x == 1)
          })]
        }
        dado <- as.data.frame(do.call("cbind", dddd))
        head(dado)
        doo <- apply(dado, 2, function(x) {
          s1 <- length(x)
          s2 <- length(unique(x))
          if (s1 == s2) {
            y <- FALSE
          }
          else {
            y <- TRUE
          }
        })
        goood <- which(doo)
        baaad <- which(!doo)
        doo2 <- names(dado)[which(doo)]
        a <- paste("y~", paste(paste("(1|", doo2, ")", 
                                     sep = ""), collapse = "+"))
        lmermodel <- lmer(as.formula(a), data = dado)
        vc <- as.data.frame(VarCorr(lmermodel))
        rownames(vc) <- vc$grp
        if (length(baaad) > 0) {
          bad.mat <- matrix(0, ncol = 5, nrow = length(baaad))
          rownames(bad.mat) <- names(baaad)
          colnames(bad.mat) <- colnames(vc)
          vc2 <- rbind(vc, bad.mat)
        }
        else {
          vc2 <- vc
        }
        were.squares <- unlist(lapply(Z, function(x) {
          if (!is.null(attributes(x$K)$diagon)) {
            y <- TRUE
          }
          else {
            y <- FALSE
          }
          return(y)
        }))
        conditionK <- length(which(were.squares)) == 
          length(Z)
        if (conditionK & (length(baaad) == 0)) {
          forced <- vc2[c(colnames(dado), "Residual"), 
                        4]
        }
        else {
          init <- vc2[c(colnames(dado), "Residual"), 
                      4]
        }
      }
    }
  }
  if (is.list(Z)) {
    if (!is.null(Z) & !is.null(W)) {
      y[which(is.na(y))] <- mean(y, na.rm = TRUE)
      misso <- which(is.na(y))
      if (length(misso) > 0) {
        y[misso] <- mean(y, na.rm = TRUE)
      }
      if (is.null(colnames(W))) {
        colnames(W) <- paste("M", 1:dim(W)[2], sep = "-")
      }
      W <- apply(W, 2, function(x) {
        vv <- which(is.na(x))
        if (length(vv) > 0) {
          mu <- mean(x, na.rm = TRUE)
          x[vv] <- mu
        }
        else {
          x <- x
        }
      })
      if (!silent) {
        cat("Estimating variance components\n")
      }
      fixed <- which(unlist(lapply(Z, function(x) {
        names(x)[1]
      })) == "X")
      random <- which(unlist(lapply(Z, function(x) {
        names(x)[1]
      })) == "Z")
      random2 <- which(names(Z) == "Z")
      if (n.PC > 0) {
        KK <- A.mat(W, shrink = FALSE)
        eig.vec <- eigen(KK)$vectors
        if (is.null(X)) {
          X <- as.matrix(rep(1, dim(KK)[1]))
        }
        Zss <- lapply(Z, function(x) {
          x[[1]]
        })
        Zssp <- as(do.call("cbind", Zss), Class = "sparseMatrix")
        X <- make.full(cbind(X, Zssp %*% (1 + as.matrix(eig.vec[, 
                                                                1:n.PC]))))
      }
      if (length(random) > 1 & method == "EMMA") {
        stop
        cat("\nError. The EMMA and GEMMA methods were design to only deal with a single variance component other than error, please select method='AI', method='NR', or method='EM' which can estimate more than one variance component\n\n")
      }
      if ((length(random) == 1 | length(random2) == 1) & 
          method == "EMMA") {
        if (length(random) > 0) {
          res <- EMMA(y = y, X = X, Z = Z[[random]][[1]], 
                      K = Z[[random]][[2]], REML = REML, silent = silent, 
                      EIGEND = EIGEND)
          names(res$u.hat) <- names(Z)
          if (!is.null(names(Z))) {
            rownames(res$var.comp) <- c(paste("Var(", 
                                              names(Z), ")", sep = ""), "Var(Error)")
          }
        }
        else {
          res <- EMMA(y = y, X = X, Z = Z[[1]], K = Z[[2]], 
                      REML = REML, silent = silent, EIGEND = EIGEND)
          names(res$u.hat) <- names(Z)
          if (!is.null(names(Z))) {
            rownames(res$var.comp) <- c(paste("Var(", 
                                              names(Z), ")", sep = ""), "Var(Error)")
          }
        }
      }
      if (method == "EM") {
        res <- EM(y = y, X = X, ETA = Z, R = R, init = init, 
                  iters = iters, REML = REML, draw = draw, silent = silent, 
                  forced = forced)
      }
      if (method == "AI") {
        if (length(Z) == 1) {
          dias <- unlist(lapply(Z[[1]], function(x) {
            if (dim(x)[1] == dim(x)[2]) {
              y <- is.diagonal.matrix(x)
            }
            else {
              y <- FALSE
            }
            return(y)
          }))
          if (length(which(dias)) == 2) {
            res <- EMMA(y = y, X = X, Z = Z[[random]][[1]], 
                        K = Z[[random]][[2]], REML = REML, silent = silent, 
                        EIGEND = EIGEND)
            names(res$u.hat) <- names(Z)
            if (!is.null(names(Z))) {
              rownames(res$var.comp) <- c(paste("Var(", 
                                                names(Z), ")", sep = ""), "Var(Error)")
            }
          }
          else {
            res <- AI(y = y, X = X, ZETA = Z, R = R, 
                      REML = REML, draw = draw, silent = silent, 
                      iters = iters, constraint = constraint, 
                      init = init, sherman = sherman, che = FALSE, 
                      EIGEND = EIGEND, Fishers = Fishers, forced = forced)
          }
        }
        else {
          res <- AI(y = y, X = X, ZETA = Z, R = R, REML = REML, 
                    draw = draw, silent = silent, iters = iters, 
                    constraint = constraint, init = init, sherman = sherman, 
                    che = FALSE, EIGEND = EIGEND, Fishers = Fishers, 
                    forced = forced)
        }
      }
      if (method == "NR") {
        if (length(Z) == 1) {
          dias <- unlist(lapply(Z[[1]], function(x) {
            if (dim(x)[1] == dim(x)[2]) {
              y <- is.diagonal.matrix(x)
            }
            else {
              y <- FALSE
            }
            return(y)
          }))
          if (length(which(dias)) == 2) {
            res <- EMMA(y = y, X = X, Z = Z[[random]][[1]], 
                        K = Z[[random]][[2]], REML = REML, silent = silent, 
                        EIGEND = EIGEND)
            names(res$u.hat) <- names(Z)
            if (!is.null(names(Z))) {
              rownames(res$var.comp) <- c(paste("Var(", 
                                                names(Z), ")", sep = ""), "Var(Error)")
            }
          }
          else {
            res <- NR(y = y, X = X, ZETA = Z, R = R, 
                      REML = REML, draw = draw, silent = silent, 
                      iters = iters, constraint = constraint, 
                      init = init, sherman = sherman, che = FALSE, 
                      Fishers = Fishers, forced = forced)
          }
        }
        else {
          res <- NR(y = y, X = X, ZETA = Z, R = R, REML = REML, 
                    draw = draw, silent = silent, iters = iters, 
                    constraint = constraint, init = init, sherman = sherman, 
                    che = FALSE, Fishers = Fishers, forced = forced)
        }
      }
      X2 <- res$X
      min.MAF = min.MAF
      n <- length(y)
      cat("\nPerforming GWAS")
      max.geno.freq = 1 - min.MAF
      W.scores <- list(NA)
      if (length(models) > 2) {
        layout(matrix(1:4, 2, 2))
      }
      else {
        layout(matrix(1:length(models), 1, length(models)))
      }
      qq <- function(scores) {
        remove <- which(scores == 0)
        if (length(remove) > 0) {
          x <- sort(scores[-remove], decreasing = TRUE)
        }
        else {
          x <- sort(scores, decreasing = TRUE)
        }
        n <- length(x)
        unif.p <- -log10(ppoints(n))
        plot(unif.p, x, pch = 16, xlab = expression(paste("Expected ", 
                                                          -log[10], "(p.value)")), ylab = expression(paste("Observed ", 
                                                                                                           -log[10], "(p.value)")), col = transp("cadetblue"), 
             main = "QQ-plot")
        lines(c(0, max(unif.p, na.rm = TRUE)), c(0, max(unif.p, 
                                                        na.rm = TRUE)), lty = 2, lwd = 2, col = "blue")
      }
      deviations <- apply(W, 2, sd)
      dev.no0 <- which(deviations > 0)
      W <- W[, dev.no0]
      for (u in 1:length(models)) {
        model <- models[u]
        cat(paste("\nRunning", model, "model"))
        ZO <- diag(dim(W)[1])
        step2 <- score.calc(marks = colnames(W), y = y, 
                            Z = ZO, X = X2, K = res$K, ZZ = res$Z, M = W, 
                            Hinv = res$V.inv, ploidy = ploidy, model = model, 
                            min.MAF = min.MAF, max.geno.freq = max.geno.freq, 
                            silent = silent, P3D = P3D)
        W.scores[[u]] <- as.matrix(step2$score)
        rownames(W.scores[[u]]) <- colnames(W)
        if (!is.null(map)) {
          dd <- W.scores[[u]]
          ffr <- fdr(dd, fdr.level = fdr.level)$fdr.10
          non.dup <- which(!duplicated(map$Locus))
          map2 <- map[non.dup, ]
          rownames(map2) <- map2$Locus
          intro <- intersect(rownames(map2), rownames(dd))
          choco <- which(colnames(map2) == "Chrom")
          if (length(intro) > 0 & length(choco) > 0) {
            map3 <- map2[intro, ]
            dd2 <- as.matrix(dd[intro, ])
            map3$p.val <- dd[intro, ]
            if (is.null(manh.col)) {
              col.scheme <- rep((transp(c("cadetblue", 
                                          "red"))), 30)
            }
            else {
              col.scheme <- rep(manh.col, 30)
            }
            layout(matrix(c(1, 2, 2), 1, 3))
            if (gwas.plots) {
              qq(step2$score)
              plot(dd2, bty = "n", col = col.scheme[factor(map3$Chrom, 
                                                           levels = unique(map3$Chrom, na.rm = TRUE))], 
                   xaxt = "n", xlab = "Chromosome", ylab = expression(paste(-log[10], 
                                                                            "(p.value)")), pch = 20, cex = 2.5, 
                   las = 2)
              init.mrks <- apply(data.frame(unique(map3$Chrom)), 
                                 1, function(x, y) {
                                   z <- which(y == x)[1]
                                   return(z)
                                 }, y = map3$Chrom)
              fin.mrks <- apply(data.frame(unique(map3$Chrom)), 
                                1, function(x, y) {
                                  z <- which(y == x)
                                  z2 <- z[length(z)]
                                  return(z2)
                                }, y = map3$Chrom)
              inter.mrks <- init.mrks + ((fin.mrks - 
                                            init.mrks)/2)
              axis(side = 1, at = inter.mrks, labels = paste("Chr", 
                                                             unique(map3$Chrom), sep = ""), cex.axis = 0.5)
              abline(h = ffr, col = "slateblue4", lty = 3, 
                     lwd = 2)
              legend("topright", legend = paste("FDR(", 
                                                fdr.level, ")=", round(ffr, 2), sep = ""), 
                     bty = "n", lty = 3, lwd = 2, col = "slateblue4", 
                     cex = 0.8)
            }
          }
          else {
            cat("\nError found! There was no markers in common between the column names of the W matrix \nand the map you provided. Please make sure that your data frame has names \n'Chrom' and 'Locus' to match correctly your map and markers tested. Plotting all markers.\n")
            map3 <- NULL
            layout(matrix(1:2, 1, 2))
            qq(step2$score)
            plot(step2$score, col = transp("cadetblue", 
                                           0.6), pch = 20, xlab = "Marker index", 
                 ylab = expression(paste(-log[10], "(p.value)")), 
                 main = paste(model, "model"), bty = "n", 
                 cex = 1.5)
          }
        }
        else if (is.null(map) & gwas.plots) {
          ffr <- fdr(step2$score, fdr.level = fdr.level)$fdr.10
          layout(matrix(1:2, 1, 2))
          qq(step2$score)
          map3 <- NULL
          plot(step2$score, col = transp("cadetblue", 
                                         0.6), pch = 20, xlab = "Marker index", ylab = expression(paste(-log[10], 
                                                                                                        "(p.value)")), main = paste(model, "model"), 
               bty = "n", cex = 1.5)
          abline(h = ffr, col = "slateblue4", lty = 3, 
                 lwd = 2)
        }
      }
      names(W.scores) <- models
      res$W.scores <- W.scores
      res$W <- W
      if (!is.null(map)) {
        res$map <- map3
      }
      res$method <- method
    }
  }
  if (is.list(Z)) {
    if ((!is.null(Z) & is.null(W))) {
      if (!silent) {
        cat("Estimating variance components\n")
      }
      fixed <- which(unlist(lapply(Z, function(x) {
        names(x)[1]
      })) == "X")
      random <- which(unlist(lapply(Z, function(x) {
        names(x)[1]
      })) == "Z")
      random2 <- which(names(Z) == "Z")
      if (length(random) > 1 & method == "EMMA") {
        stop
        cat("\nError, The EMMA and GEMMA methods were design to only deal with a single variance component other than error, please select method='AI', method='NR', or method='EM' which can estimate more than one variance component\n\n")
      }
      if ((length(random) == 1 | length(random2) == 1) & 
          method == "EMMA") {
        if (length(random) > 0) {
          res <- EMMA(y = y, X = X, Z = Z[[random]][[1]], 
                      K = Z[[random]][[2]], REML = REML, silent = silent, 
                      EIGEND = EIGEND)
          names(res$u.hat) <- names(Z)
          if (!is.null(names(Z))) {
            rownames(res$var.comp) <- c(paste("Var(", 
                                              names(Z), ")", sep = ""), "Var(Error)")
          }
        }
        else {
          res <- EMMA(y = y, X = X, Z = Z[[1]], K = Z[[2]], 
                      REML = REML, silent = silent, EIGEND = EIGEND)
          names(res$u.hat) <- names(Z)
          if (!is.null(names(Z))) {
            rownames(res$var.comp) <- c(paste("Var(", 
                                              names(Z), ")", sep = ""), "Var(Error)")
          }
        }
      }
      if (method == "EM") {
        res <- EM(y = y, X = X, ETA = Z, R = R, init = init, 
                  iters = iters, REML = REML, draw = draw, silent = silent, 
                  forced = forced)
      }
      if (method == "AI") {
        if (length(Z) == 1) {
          dias <- unlist(lapply(Z[[1]], function(x) {
            if (dim(x)[1] == dim(x)[2]) {
              y <- is.diagonal.matrix(x)
            }
            else {
              y <- FALSE
            }
            return(y)
          }))
          if (length(which(dias)) == 2) {
            res <- EMMA(y = y, X = X, Z = Z[[random]][[1]], 
                        K = Z[[random]][[2]], REML = REML, silent = silent, 
                        EIGEND = EIGEND)
            names(res$u.hat) <- names(Z)
            if (!is.null(names(Z))) {
              rownames(res$var.comp) <- c(paste("Var(", 
                                                names(Z), ")", sep = ""), "Var(Error)")
            }
          }
          else {
            res <- AI(y = y, X = X, ZETA = Z, R = R, 
                      REML = REML, draw = draw, silent = silent, 
                      iters = iters, constraint = constraint, 
                      init = init, sherman = sherman, che = FALSE, 
                      EIGEND = EIGEND, Fishers = Fishers, gss = gss, 
                      forced = forced)
          }
        }
        else {
          res <- AI(y = y, X = X, ZETA = Z, R = R, REML = REML, 
                    draw = draw, silent = silent, iters = iters, 
                    constraint = constraint, init = init, sherman = sherman, 
                    che = FALSE, EIGEND = EIGEND, Fishers = Fishers, 
                    gss = gss, forced = forced)
        }
      }
      if (method == "NR") {
        if (length(Z) == 1) {
          dias <- unlist(lapply(Z[[1]], function(x) {
            if (dim(x)[1] == dim(x)[2]) {
              y <- is.diagonal.matrix(x)
            }
            else {
              y <- FALSE
            }
            return(y)
          }))
          if (length(which(dias)) == 2) {
            res <- EMMA(y = y, X = X, Z = Z[[random]][[1]], 
                        K = Z[[random]][[2]], REML = REML, silent = silent, 
                        EIGEND = EIGEND)
            names(res$u.hat) <- names(Z)
          }
          else {
            res <- NR(y = y, X = X, ZETA = Z, R = R, 
                      REML = REML, draw = draw, silent = silent, 
                      iters = iters, constraint = constraint, 
                      init = init, sherman = sherman, che = FALSE, 
                      Fishers = Fishers, forced = forced)
          }
        }
        else {
          res <- NR(y = y, X = X, ZETA = Z, R = R, REML = REML, 
                    draw = draw, silent = silent, iters = iters, 
                    constraint = constraint, init = init, sherman = sherman, 
                    che = FALSE, Fishers = Fishers, forced = forced)
        }
      }
      res$method <- method
      res$maxim <- REML
      res$W <- W
    }
  }
  if ((is.null(X) & is.null(Z) & is.null(W))) {
    res <- lm(y ~ 1)
  }
  if ((!is.null(X) & is.null(Z) & is.null(W))) {
    res <- lm(y ~ X - 1)
  }
  class(res) <- c("mmer")
  layout(matrix(1, 1, 1))
  return(res)
}