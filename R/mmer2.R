mmer2 <- function(fixed, random, G=NULL, R=NULL, W=NULL, method="NR", REML=TRUE, MVM=FALSE, iters=20, draw=FALSE, init=NULL, data, family=gaussian, silent=FALSE, constraint=TRUE, sherman=FALSE, EIGEND=FALSE, gss=TRUE, forced=NULL, map=NULL, fdr.level=0.05, manh.col=NULL, min.n=FALSE, gwas.plots=TRUE, n.cores=1, tolpar = 1e-06, tolparinv = 1e-06){
  if(missing(data)){
    data <- environment(fixed)
    data2 <- environment(random)
    nodata <-TRUE
  }else{nodata=FALSE} 
  
  if(missing(random)){
    stop("Please use 'lm' for fixed effect models", call. = FALSE)
  }
  if(!is.null(G) & method == "EM"){
    cat("With var-cov structures (G) present you may want to try the AI or NR algorithm.\n\n")
  }
  ############## impute data
  data2 <- data
  for(i in 1:dim(data2)[2]){
    x <- data2[,i]
    if(is.numeric(x)){
      isNA <- which(is.na(x))
      if(length(isNA) >0){x[isNA] <- mean(x,na.rm=TRUE)}
    }else if(is.factor(x)){
      isNA <- which(is.na(x))
      if(length(isNA) >0){x[isNA] <- as.factor(names(which(table(x) == max(table(x)))[1]))}
    }else if(is.character(x)){
      isNA <- which(is.na(x))
      if(length(isNA) >0){x[isNA] <- as.character(names(which(table(x) == max(table(x)))[1]))}
    }
    data2[,i] <- x
  }
  ######################
  
  mf <- try(model.frame(fixed, data = data2, na.action = na.pass), silent = TRUE)
  mfna <- try(model.frame(fixed, data = data, na.action = na.pass), silent = TRUE)
  if (class(mf) == "try-error") {
    stop("Please provide the 'data' argument for your specified variables", call. = FALSE)
  }
  mf <- eval(mf, parent.frame())
  mfna <- eval(mfna, parent.frame())
  
  #which(!duplicated(t(mfna)))
  # response Y
  yvar <- model.response(mfna)
  
  #yvar <- gsub(" ", "", as.character(fixed[2]))
  ### Xb in 'formula'
  X <- model.matrix(fixed, mf)

  #xvar <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[+]")[[1]])
  ### Zu in formula
  
  if(!is.null(random)){
    if(nodata){
      V <- try(model.frame(random, data = data2, na.action = na.pass), silent = TRUE)
      if (class(V) == "try-error") {
        stop("Please provide the 'data' argument for your specified variables", call. = FALSE)
      }#V <- model.frame(random, data = data2, na.action = na.pass)
      V <- eval(V, parent.frame())
    }else{
      V <- try(model.frame(random, data = data, na.action = na.pass), silent = TRUE)
      if (class(V) == "try-error") {
        stop("Please provide the 'data' argument for your specified variables", call. = FALSE)
      }#V <- model.frame(random, data = data, na.action = na.pass)
      V <- eval(V, parent.frame()) 
    }
    
    zvar.names <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+]")[[1]])
    #zvar.names <- zvar.names[which(!duplicated(zvar.names))]
    #print(zvar.names)
    zvar <- V #names(V)
    for(i in 1:dim(zvar)[2]){
      zvar[,i] <- as.factor(zvar[,i])
    }
    #zvar <- apply(zvar,2,as.factor)
    #zvar.names <- names(V)
    #zvar <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+]")[[1]])
    #varsss <- c(xvar,zvar)
    Z <- list()
    for(i in 1:length(zvar.names)){
      ## incidence matrix
      vara <- zvar.names[i]
      # data.frame(factor(V[,vara],levels=V[,vara],ordered=T))
      zi <- model.matrix(as.formula(paste("~",vara,"-1")),zvar)

      ## var-cov matrix
      ww <- which(names(G) %in% vara)
      if(length(ww) > 0){# K was provided
        ## just if there's a K matrix we make sure to be using the real names and no the model.matrix ones
        colnames(zi) <- levels(V[,vara])
        #########
        uuuz <- levels(as.factor(colnames(zi))) # order of Z
        uuuk <- attr(G[[ww]],"dimnames")[[1]] # order of K
        inte <- intersect(uuuz,uuuk)
        if(length(inte)==length(uuuz)){ # the names were the same in Z and K
          ki <- G[[ww]][uuuz,uuuz]
        }else{ # no intersection between z and k names
          cat(paste("\nNames of Z and K for random effect",vara,"are not the same. \nMake sure they are in the correct order."))
          ki <- G[[ww]] 
        }
        
      }else{ # was not provided, we create a diagonal
        ki <- diag(dim(zi)[2])
      }
      elem <- list(Z=zi, K=ki)
      Z[[i]] <- elem
    }
    names(Z) <- zvar.names
    
    res <- mmer(Y=yvar, X=X, Z=Z, R=R, W=W, method=method, REML=REML, iters=iters, draw=draw, init=init, silent=silent, constraint=constraint, sherman=sherman, EIGEND=EIGEND, gss=gss, forced=forced, map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=gwas.plots,n.cores=n.cores, MVM=MVM,tolpar = tolpar, tolparinv = tolparinv)
    
  }else{###only fixed effects
    res <- glm(yvars~X, family=family)
  }
  #########
  return(res)
}
