EM <- function (y, X = NULL, ETA = NULL, R = NULL, init = NULL, iters = 50, 
          REML = TRUE, draw = TRUE, silent = FALSE) 
{
  if (is.null(X) & is.null(ETA)) {
    tn = length(y)
    xm <- matrix(1, tn, 1)
    output <- lm(y ~ xm - 1)
  }
  else {
    if (is.list(ETA)) {
      if (is.list(ETA[[1]])) {
        ETA = ETA
      }
      else {
        ETA = list(ETA)
      }
    }
    else {
      stop
      cat("\nThe random effects need to be provided in a list format, please see examples")
    }
    if (is.null(X) & !is.null(ETA)) {
      ETA0 <- list(X = matrix(rep(1, times = length(y)), 
                              ncol = 1), K = diag(dim(matrix(rep(1, times = length(y))))[2]))
      ETA <- c(list(ETA0), ETA)
    }
    if (!is.null(X) & !is.null(ETA)) {
      if (is.list(X)) {
        if (is.list(X[[1]])) {
          ETA0 <- list(X = X[[1]][[1]], K = diag(dim(X[[1]][[1]])[2]))
          ETA <- c(list(ETA0), ETA)
        }
        else {
          ETA0 <- list(X = X[[1]], K = diag(dim(X[[1]])[2]))
          ETA <- c(list(ETA0), ETA)
        }
      }
      else {
        X = X
        K <- diag(dim(X)[2])
        ETA0 <- list(X = X, K = K)
        ETA <- c(list(ETA0), ETA)
      }
    }
    if (is.null(R)) {
      R <- diag(length(y))
    }
    fixed <- which(unlist(lapply(ETA, function(x) {
      names(x)[1]
    })) == "X")
    if (length(fixed) == 0) {
      X <- matrix(rep(1, times = length(y)), ncol = 1)
      K <- diag(dim(X)[2])
      ETA0 <- list(X = X, K = K)
      ETA <- c(list(ETA0), ETA)
    }
    ETA <- lapply(ETA, function(x) {
      if (length(x) == 1) {
        x[[2]] <- diag(dim(x[[1]])[2])
      }
      else {
        x <- x
      }
      return(x)
    })
    eta.or <- ETA
    eta.or <- lapply(eta.or, function(x) {
      lapply(x, as.matrix)
    })
    ETA2 <- ETA
    y2 <- y
    good <- which(!is.na(y))
    ETA <- lapply(ETA, function(x, good) {
      x[[1]] <- x[[1]][good, ]
      x[[2]] <- x[[2]]
      return(x)
    }, good = good)
    y <- y[good]
    R <- R[good, good]
    nran <- length(which(unlist(lapply(ETA, function(x) {
      names(x)[1]
    })) == "Z"))
    var.com <- vector(mode = "list", length = nran + 1)
    if (is.null(init)) {
      var.y = var(y, na.rm = TRUE)
      var.com <- lapply(var.com, function(x) {
        x = var.y/nran
      })
    }
    var.com.sto <- var.com
    var.com[[length(var.com)]] = var.y
    ETA <- lapply(ETA, function(x) {
      lapply(x, as.matrix)
    })
    if (!silent) {
      count <- 0
      tot <- iters
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    varE = var(y, na.rm = TRUE)
    conv = 0
    wi = 0
    change = 1
    while (conv == 0) {
      wi = wi + 1
      if (!silent) {
        count <- count + 1
      }
      axs = lapply(var.com, function(x) {
        varE/x
      })
      CC = numeric()
      for (j in 1:length(ETA)) {
        prov <- numeric()
        for (k in 1:length(ETA)) {
          if (j == k & names(ETA[[j]])[1] != "X") {
            res <- crossprod(ETA[[j]][[1]], ETA[[k]][[1]]) + 
              (as.vector(axs[[j - 1]]) * solve(as.matrix(ETA[[j]][[2]])))
          }
          else {
            res <- crossprod(ETA[[j]][[1]], ETA[[k]][[1]])
          }
          prov <- cbind(prov, res)
        }
        CC <- rbind(CC, prov)
      }
      l <- lapply(ETA, function(x, y) {
        t(x[[1]]) %*% y
      }, y = y)
      l2 <- as.matrix(unlist(l))
      CInv <- solve(CC)
      thetaHat <- CInv %*% l2
      nn <- lapply(ETA, function(x) {
        dim(x[[1]])[2]
      })
      nn <- lapply(nn, function(x) {
        if (is.null(x)) {
          x = 1
        }
        else {
          x = x
        }
        return(x)
      })
      pairs.a = list(NA)
      for (h in 1:length(nn)) {
        pairs.a[[h]] <- ((sum(unlist(nn[1:h])) - (unlist(nn[h]) - 
                                                    1)):sum(unlist(nn[1:h])))
      }
      now <- 0
      for (f in 1:length(pairs.a)) {
        now <- now + crossprod(as.matrix(ETA[[f]][[1]]) %*% 
                                 thetaHat[pairs.a[[f]], ], y)
      }
      varE = (crossprod(y) - now)/(length(y) - nn[[1]])
      indexK <- which(unlist(lapply(ETA, function(x) {
        names(x)[1]
      })) == "Z")
      track <- var.com
      for (k in 1:(length(var.com) - 1)) {
        rrr <- indexK[k]
        Kinv <- solve(ETA[[rrr]][[2]])
        var.com[[k]] = ((t(thetaHat[pairs.a[[rrr]], ]) %*% 
                           Kinv %*% thetaHat[pairs.a[[rrr]], ]) + matrixcalc::matrix.trace(Kinv %*% 
                                                                                             CInv[pairs.a[[rrr]], pairs.a[[rrr]]] * as.numeric(varE)))/nrow(ETA[[rrr]][[2]])
      }
      var.com[[length(var.com)]] = varE
      for (k in 1:length(var.com.sto)) {
        var.com.sto[[k]] <- c(var.com.sto[[k]], var.com[[k]])
      }
      lege2 <- list()
      for (k in 1:length(var.com)) {
        if (k == length(var.com)) {
          lege2[[k]] <- paste("Var(e):")
        }
        else {
          lege2[[k]] <- paste("Var(u", k, "):", sep = "")
        }
      }
      if (draw) {
        ylim <- max(unlist(var.com.sto), na.rm = TRUE)
        my.palette <- RColorBrewer::brewer.pal(7, "Accent")
        layout(matrix(1, 1, 1))
        plot(var.com.sto[[1]], ylim = c(0, ylim), type = "l", 
             col = my.palette[1], lwd = 3, xlim = c(0, iters), 
             xaxt = "n", las = 2, main = "Expectation-Maximization algorithm results", 
             ylab = "Value of the variance component", xlab = "Iteration to be processed to reach convergence", 
             cex.axis = 0.8)
        axis(1, las = 1, at = 0:10000, labels = 0:10000, 
             cex.axis = 0.8)
        for (t in 1:length(var.com)) {
          lines(var.com.sto[[t]], col = my.palette[t], 
                lwd = 3)
        }
        ww <- length(var.com.sto[[1]])
        lege <- list()
        for (k in 1:length(var.com)) {
          if (k == length(var.com)) {
            lege[[k]] <- paste("Var(e):", round(var.com.sto[[k]][ww], 
                                                4), sep = "")
          }
          else {
            lege[[k]] <- paste("Var(u", k, "):", round(var.com.sto[[k]][ww], 
                                                       4), sep = "")
          }
        }
        legend("topright", bty = "n", col = my.palette, 
               lty = 1, legend = unlist(lege), cex = 0.75)
      }
      if (!silent) {
        setTxtProgressBar(pb, (count/tot))
      }
      if (wi > 1) {
        change = abs(sum(unlist(var.com) - unlist(track)))
      }
      if (change < 1e-05 | wi == iters) {
        conv = 1
        if (!silent) {
          setTxtProgressBar(pb, (tot/tot))
        }
        if (wi == iters) {
          cat("\nMaximum number of iterations reached with no convergence using the EM algorithm, look at the variance components change over iterations (plot) and be cautious using the variance components estimated if they don't look steady")
        }
      }
    }
    xm <- as.matrix(ETA[[1]][[1]])
    zvar <- which(unlist(lapply(ETA, function(x) {
      names(x)[1]
    })) == "Z")
    V <- matrix(0, nrow = dim(xm)[1], ncol = dim(xm)[1])
    V2 <- R * (as.numeric(varE))
    for (k in zvar) {
      V2 <- V2 + tcrossprod(ETA[[k]][[1]], ETA[[k]][[1]] %*% 
                              (ETA[[k]][[2]] * as.numeric(var.com[k - 1])))
      V <- V + tcrossprod(ETA[[k]][[1]], ETA[[k]][[1]] %*% 
                            (ETA[[k]][[2]])) + (R * (as.numeric(varE)/as.numeric(var.com[k - 
                                                                                           1])))
    }
    Vinv <- solve(V)
    Vinv2 <- solve(V2)
    X <- ETA[[fixed]][[1]]
    xvx = t(X) %*% Vinv %*% X
    xvxi = solve(xvx)
    pm = Vinv - Vinv %*% X %*% xvxi %*% t(X) %*% Vinv
    if (REML == TRUE) {
      logL = -0.5 * (log(det(V)) + log(det(xvx)) + t(y) %*% 
                       pm %*% y)
    }
    else {
      logL = -0.5 * (log(det(V)) + t(y) %*% pm %*% y)
    }
    xvar <- which(unlist(lapply(ETA, function(x) {
      names(x)[1]
    })) == "X")
    beta <- solve(crossprod(ETA[[xvar]][[1]], Vinv %*% ETA[[xvar]][[1]])) %*% 
      crossprod(ETA[[xvar]][[1]], Vinv %*% y)
    u <- list()
    for (k in zvar) {
      u[[k]] <- ETA[[k]][[2]] %*% t(ETA[[k]][[1]]) %*% 
        Vinv %*% (y - ETA[[xvar]][[1]] %*% beta)
    }
    u <- u[zvar]
    Var.u.hat <- list()
    PEV.u.hat <- list()
    for (i in zvar) {
      Var.u.hat[[i]] <- (as.numeric(var.com[i - 1])^2) * 
        (crossprod(ETA[[i]][[1]] %*% ETA[[i]][[2]], pm) %*% 
           (ETA[[i]][[1]] %*% ETA[[i]][[2]]))
      PEV.u.hat[[i]] <- as.numeric(var.com[i - 1]) * ETA[[i]][[2]] - 
        Var.u.hat[[i]]
    }
    Var.u.hat <- Var.u.hat[zvar]
    PEV.u.hat <- PEV.u.hat[zvar]
    Var.beta.hat <- solve(crossprod(xm, Vinv2 %*% xm))
    AIC = (-2 * logL) + (2 * dim(X)[2])
    BIC = (-2 * logL) + (log(length(y)) * dim(X)[2])
    fitted.y <- eta.or[[1]][[1]] %*% beta
    fitted.u <- 0
    for (h in 2:length(eta.or)) {
      fitted.y <- fitted.y + (eta.or[[h]][[1]] %*% u[[h - 
                                                        1]])
      if (h > 1) {
        fitted.u <- fitted.u + (eta.or[[h]][[1]] %*% 
                                  u[[h - 1]])
      }
    }
    fitted.y.good <- fitted.y[good]
    residuals2 <- y - (ETA[[1]][[1]] %*% beta)
    residuals3 <- y - fitted.y[good]
    out1 <- unlist(var.com)
    names(out1) <- unlist(lege2)
    rownames(beta) <- colnames(ETA[[1]][[1]])
    for (i in 2:length(ETA)) {
      rownames(u[[i - 1]]) <- colnames(ETA[[i]][[1]])
    }
    output <- list(var.comp = out1, V.inv = Vinv, u.hat = u, 
                   Var.u.hat = Var.u.hat, PEV.u.hat = PEV.u.hat, beta.hat = beta, 
                   Var.beta.hat = Var.beta.hat, LL = logL, AIC = AIC, 
                   BIC = BIC, X = ETA[[xvar]][[1]], fitted.y = fitted.y, 
                   fitted.u = fitted.u, residuals = residuals2, cond.residuals = residuals3, 
                   fitted.y.good = fitted.y.good)
  }
  return(output)
}