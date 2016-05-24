mmer2 <- function(fixed=NULL, random=NULL, G=NULL, R=NULL, W=NULL, method="NR", REML=TRUE, iters=50, draw=FALSE, init=NULL, data=NULL, family=gaussian, silent=FALSE, constraint=TRUE, sherman=FALSE, MTG2=FALSE, gss=FALSE, forced=NULL, map=NULL, fdr.level=0.05, manh.col=NULL, min.n=TRUE, gwas.plots=TRUE){
  if(!is.null(G) & method == "EM"){
    cat("With var-cov structures (G) present you may want to try the AI algorithm.\n\n")
  }
  ### Response "y"
  yvar <- gsub(" ", "", as.character(fixed[2]))
  ### Xb in 'formula'
  xvar <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[+]")[[1]])
  ### Zu in formula
  if(!is.null(random)){
    zvar <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+]")[[1]])
    varsss <- c(xvar,zvar)
  }else{varsss <- xvar}
  varsss <- varsss[which(varsss != "1")]
  doto <- grep(":",varsss)
  if(length(doto) >0){varsss <- varsss[-doto]}
  ### FIX THE DATA IF NA's EXIST
  good <- list()
  for(i in 1:length(varsss)){
    xxx <- data[,varsss[i]]
    good[[i]] <- which(!is.na(xxx))
  }
  
  if(length(varsss) == 1){
    keep <- as.vector(unlist(good))
  }else{
    keep <- good[[1]]
    for(j in 1:length(good)){
      keep <- intersect(keep,good[[j]])
    }
  }
  data <- data[keep,] # only good data in model
  ##############################
  if(length(xvar %in% "1") == 1){
    X <- as.matrix(model.matrix(as.formula(paste("~ ", paste(c(xvar), collapse="+"))), data=data))
  }else{
    X <- as.matrix(model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=data))
  }
  
  ### random and G
  if(!is.null(random)){ # ============== MIXED MODEL =====================
    #zvar <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+]")[[1]])
    # generalized mixed model, carefull, NA's dissapear
    #mox <- glm(data[,yvar] ~ 1, family = family)
    yvars <- data[,yvar] #mox$family$linkfun(mox$y)
    #
    if(is.null(random)){
      Z=NULL
    }else{
      Z <- list()
      for(i in 1:length(zvar)){ # do it factor
        ## incidence matrix
        vara <- zvar[i]
        data2 <- data.frame(apply(data,2,as.factor))
        zi <- model.matrix((as.formula(paste("~ -1 + ", vara))), data=data2)
        colnames(zi) <- gsub(vara,"",colnames(zi))
        if(dim(zi)[2] == length(yvars) & min.n == TRUE){ # lmer error
          stop("Error: number of levels of each grouping factor must be < number of observations")
        }else{ # right way to specify the random effects, keep going
          ## var-cov matrix
          ww <- which(names(G) %in% vara)
          if(length(ww) > 0){# was provided
            ki <- G[[ww]]
          }else{ # was not provided, we create a diagonal
            ki <- diag(dim(zi)[2])
          }
          elem <- list(Z=zi, K=ki)
          Z[[i]] <- elem
        }
      }
      names(Z) <- zvar
      ### fit the model using the real function mmer2 
      res <- mmer(y=yvars, X=X, Z=Z, R=R, W=W, method=method, REML=REML, iters=iters, draw=draw, init=init, silent=silent, constraint=constraint, sherman=sherman, MTG2=MTG2, gss=gss, forced=forced, map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=gwas.plots)
      #rownames(res$var.comp) <- c(zvar,"Error")
    }
  }else{ # ================== JUST FIXED =======================
    res <- glm(yvars~X, family=family)
  }
  class(res)<-c("mmer")
  return(res)
}
