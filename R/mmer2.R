mmer2 <- function (fixed = NULL, random = NULL, G = NULL, R = NULL, method = "AI", 
                   REML = TRUE, iters = 50, draw = FALSE, init = NULL, data = NULL, 
                   family = gaussian, silent = FALSE) 
{
  yvar <- gsub(" ", "", as.character(fixed[2]))
  xvar <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[+]")[[1]])
  if (xvar %in% "1") {
    X <- as.matrix(model.matrix(as.formula(paste("~ ", paste(c(xvar), 
                                                             collapse = "+"))), data = data))
  }
  else {
    X <- as.matrix(model.matrix(as.formula(paste("~ ", paste(c("1", 
                                                               xvar), collapse = "+"))), data = data))
  }
  if (!is.null(random)) {
    zvar <- gsub(" ", "", strsplit(as.character(random[2]), 
                                   split = "[+]")[[1]])
    yvars <- data[, yvar]
    if (is.null(random)) {
      Z = NULL
    }
    else {
      Z <- list()
      for (i in 1:length(zvar)) {
        vara <- zvar[i]
        data2 <- data.frame(apply(data, 2, as.factor))
        zi <- model.matrix((as.formula(paste("~ -1 + ", 
                                             vara))), data = data2)
        if (dim(zi)[2] == length(yvars)) {
          stop("Error: number of levels of each grouping factor must be < number of observations")
        }
        else {
          ww <- which(names(G) %in% vara)
          if (length(ww) > 0) {
            ki <- G[[ww]]
          }
          else {
            ki <- diag(dim(zi)[2])
          }
          elem <- list(Z = zi, K = ki)
          Z[[i]] <- elem
        }
      }
      res <- mmer(y = yvars, X = X, Z = Z, R = R, method = method, 
                  REML = REML, iters = iters, draw = draw, init = init, 
                  silent = silent)
    }
  }
  else {
    res <- glm(yvars ~ X, family = family)
  }
  class(res) <- c("mmer")
  return(res)
}