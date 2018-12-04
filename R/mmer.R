mmer <- function(fixed, random, rcov, data, weights, 
                 iters=20, tolpar = 1e-03, tolparinv = 1e-06, 
                 init=NULL, constraints=NULL, method="NR", 
                 getPEV=TRUE,
                 na.method.X="exclude",
                 na.method.Y="exclude",
                 return.param=FALSE, 
                 date.warning=TRUE,
                 verbose=TRUE,reshape.output=TRUE){
  
  my.year <- 2019
  my.month <- 3 #month when the user will start to get notifications the 1st day of next month
  ### if my month = 3, user will start to get notification in april 1st (next month)
  datee <- Sys.Date()
  year.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
  your.year <- year.mo.day[1]
  your.month <- year.mo.day[2]
  ## if your month is greater than my month you are outdated
  if(date.warning){
    if(your.month > my.month & your.year >= my.year){
      # error if your month is greater and your year is smaller
      cat("Version out of date. Please update sommer to the newest version using:\ninstall.packages('sommer') in a new session\n Use the 'date.warning' argument to disable the warning message.")
    }
  }
  
  if(missing(data)){
    data <- environment(fixed)
    if(!missing(random)){
      data2 <- environment(random)
    }
    nodata <-TRUE
    cat("data argument not provided \n")
  }else{nodata=FALSE} 
  
  if(missing(rcov)){
    rcov = as.formula("~units")
  }
  #################
  ## do the needed for na.method.Y and na.method.X
  dataor <- data
  provdat <- subdata(data, fixed=fixed, na.method.Y = na.method.Y,na.method.X=na.method.X)
  # print(str(provdat))
  data <- provdat$datar
  # randomization <- sample(1:nrow(data))
  # data <- data[randomization,]
  #################
  data$units <- levels(as.factor(paste("u",1:nrow(data),sep="")))
  #################
  ## get Y matrix
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (class(mfna) == "try-error") {
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mfna <- eval(mfna, parent.frame())
  yvar <- as.matrix(model.response(mfna))
  nt <- ncol(yvar)
  if(nt==1){colnames(yvar) <- response}
  
  #################
  ## get Zs and Ks
  if(!missing(random)){
  yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]]
  rtermss <- apply(data.frame(yuyu),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  zs <- list()
  ks <- list()
  ges <- list()
  gesI <- list()
  re_namel1 <- list()
  # print(str(data))
  # print(str(model.matrix(~id-1,data)))
  # print(str(as(model.matrix(~id-1,data), Class = "sparseMatrix")))
  for(u in 1:length(rtermss)){
    checkvs <- grep("vs\\(",rtermss[u])
    check2d <- grep("spl2D\\(",rtermss[u])
    if(length(checkvs)>0){ ## if this term is a variance structure
      # print(mm)
      ff <- eval(parse(text = rtermss[u]),envir = data) 
      # print(nrow(ff$Z[[1]]))
      # print(length(provdat$good))
      if(nrow(ff$Z[[1]]) != length(provdat$good)){
        ## if the incidence matrix is different than the size of good very likely the user provided a matrix
        ff$Z <- lapply(ff$Z,function(xxx){as.matrix(xxx[provdat$good,])})
      }
      # print(rtermss[u])
      # print(ff$Gtc)
      re_namel1[[u]] <- ff$re_name
      zs[[u]] <- ff$Z
      ks[[u]] <- ff$K
      if(is.null(ff$Gt)){ ## initial vc if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          if(ff$typevc[k] == 2){div=2}else{div=1}## divisor
          mml[[k]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt))/div 
        }
        ges[[u]] <- mml
      }else{ges[[u]] <- rep(list(ff$Gt),length(ff$Z))}
      # print(str(ff))
      # print(ff$Gtc)
      if(is.null(ff$Gtc)){ ## contraints if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          mm <- matrix(as.vector(ff$typevc[k]),nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
          mml[[k]] <- mm
        }
        gesI[[u]] <- mml
      }else{
        ff$Gtc[lower.tri(ff$Gtc)] <- 0
        gesI[[u]] <- rep(list(ff$Gtc),length(ff$Z))
      }
      # print(ff$Gtc)
      # print(gesI[[u]])
    }else{ ## if is a normal term
      if(length(check2d) > 0){
        ff <- eval(parse(text = rtermss[u]),envir = data)  
        re_namel1[[u]] <- names(ff)
        zs[[u]] <- ff
        ks[[u]] <- lapply(ff,function(x){pki <- diag(ncol(x));colnames(pki)<-rownames(pki) <- colnames(x);return(pki)})
        gepki <- ( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)
        ges[[u]] <- rep(list(gepki),length(ff))
        mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
        gesI[[u]] <- rep(list(mm),length(ff))
      }else{
        zpp <- model.matrix(as.formula(paste("~",rtermss[u],"-1")), data=data)
        zs[[u]] <- list(zpp)
        nu <- ncol(zpp)
        ks[[u]] <- list(diag(1, nu,nu))
        re_namel1[[u]] <- list(rtermss[u])
        ges[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)) 
        mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
        gesI[[u]] <- list(mm)
      }
    }
  }
  Z <- unlist(zs,recursive=FALSE)
  K <- unlist(ks,recursive=FALSE)
  # print(ges)
  ges <- unlist(ges,recursive=FALSE)
  gesI <- unlist(gesI,recursive=FALSE)
  re_namel1 <- unlist(re_namel1,recursive=FALSE)
  }else{re_namel1 <- character()}
  #################
  ## get Rs
  yuyur <- strsplit(as.character(rcov[2]), split = "[+]")[[1]]
  rcovtermss <- apply(data.frame(yuyur),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  rs <- list()
  gesr <- list()
  gesIr <- list()
  re_namel2 <- list()
  for(u in 1:length(rcovtermss)){
    checkvs <- grep("vs\\(",rcovtermss[u])
    if(length(checkvs)>0){ ## if this term is a variance structure
      ff <- eval(parse(text = rcovtermss[u]),envir = data)  
      # if(nrow(ff$Z) != length(provdat$good)){
        ## if the incidence matrix is different than the size of good very likely the user provided a matrix
        # ff$Z <- ff$Z[provdat$good,provdat$good]
      # }
      rs[[u]] <- ff$Z
      re_namel2[[u]] <- ff$re_name
      if(is.null(ff$Gt)){ ## initial vc if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          if(ff$typevc[k] == 2){div=2}else{div=1}## divisor
          mml[[k]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.04977728 + diag(0.02488864, nt,nt))/div 
        }
        gesr[[u]] <-  mml
      }else{gesr[[u]] <- rep(list(ff$Gt),length(ff$Z))}
      if(is.null(ff$Gtc)){ ## contraints if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          mm <- matrix(as.vector(ff$typevc[k]),nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
          mml[[k]] <- mm
        }
        gesIr[[u]] <- mml
      }else{
        ff$Gtc[lower.tri(ff$Gtc)] <- 0
        gesIr[[u]] <- rep(list(ff$Gtc),length(ff$Z))
      }
    }else{ ## if is a normal term
      rpp <- model.matrix(as.formula(paste("~",rcovtermss[u],"-1")), data=data)
      rs[[u]] <- list(rpp)
      re_namel2[[u]] <- rcovtermss[u]
      gesr[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.04977728 + diag(0.02488864, nt,nt))
      mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
      gesIr[[u]] <- list(mm)
    }
  }
  R <- unlist(rs,recursive=FALSE)
  gesr <- unlist(gesr,recursive=FALSE)
  gesIr <- unlist(gesIr,recursive=FALSE)
  re_namel2 <- unlist(re_namel2,recursive=FALSE)
  
  if(!missing(random)){
    GES <- c(ges,gesr)
    GESI <- c(gesI,gesIr)
  }else{
    GES <- c(gesr)
    GESI <- c(gesIr)
  }

  #################
  ## get Xs
  
  yuyuf <- strsplit(as.character(fixed[3]), split = "[+-]")[[1]]
  fixedtermss <- apply(data.frame(yuyuf),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  fixedtermss <- fixedtermss[which(fixedtermss != "-1")]
  if(length(fixedtermss)==1 & fixedtermss[1] == "1"){
    
  }else{
    fixedtermss <- fixedtermss[which(fixedtermss != "1")]
  }
  
  xs <- list()
  gesf <- list()
  gesIf <- list()
  
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  purefixedtermss <- apply(data.frame(fixedtermss),1,function(y){
    newy <- expi(y)
    if(length(newy) > 0){
      newy <- gsub(",.*","",newy)
    }else{newy <- y}
    return(newy)
  })
  
  if(length(fixedtermss)==1 & fixedtermss[1] == "1"){
    types <- character()
    fixedtermss <- fixedtermss[which(fixedtermss != "1")]
    purefixedtermss <- purefixedtermss[which(purefixedtermss != "1")]
  }else{
    types <- unlist(lapply(as.data.frame(data[,purefixedtermss]),class))
    names(types) <- purefixedtermss
  }
  
  # print(myformula(fixed))
  # newfixed <- as.formula(myformula(fixed))
  newfixed <- as.formula(myformula(fixed))
  mf <- try(model.frame(newfixed, data = data, na.action = na.pass), silent = TRUE)
  mf <- eval(mf, parent.frame())
  baseX <- model.matrix(newfixed, mf)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if("(Intercept)" %in% colnames(baseX) & "(Intercept)" %!in% purefixedtermss){
    fixedtermss <- c("(Intercept)",fixedtermss)
    purefixedtermss <- c("(Intercept)",purefixedtermss)
    types <- c("(Intercept)",types)
  }
  
  
  # print(head(baseX))
  
  for(u in 1:length(fixedtermss)){
    checkvs <- grep("vs\\(",fixedtermss[u])
    if(length(checkvs)>0){ ## if this term is a variance structure
      # orterm <- gsub(",.*","",expi2(fixedtermss[u]))
      # print(fixedtermss[u])
      ffx <- eval(parse(text = fixedtermss[u]),envir = data)
      # print(ffx$Gtc)
      # print(str(ffx))
       # print(head(ff))
      # xs[[u]] <- list(
      #   as.matrix(baseX[,which(gsub("(Intercept)","",colnames(baseX)) %in% colnames(ff$Z[[1]]))])
      # )
      if(purefixedtermss[u] == "(Intercept)"){
        xpp <- as.matrix(baseX[,"(Intercept)"])
        colnames(xpp) <- "(Intercept)"
      }else{
        if(types[purefixedtermss[u]] == "numeric"){ #if is numeric just look for the name of the column
          findlevs <- purefixedtermss[u]
          xpp <- as.matrix(baseX[,findlevs]); colnames(xpp) <- purefixedtermss[u]
        }else{
          findlevs <- paste0( purefixedtermss[u] , unique(as.character(na.omit(data[,purefixedtermss[u]]))) )
          xpp <- as.matrix(baseX[, which(colnames(baseX) %in% findlevs)])
        }
      }
      # print(ff$Gtc)
      xs[[u]] <- list(xpp)
      if(is.null(ffx$Gt)){ ## initial vc if user don't provide them
        gesf[[u]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)) 
      }else{gesf[[u]] <- rep(list(ffx$Gt),length(ffx$Z))}
      if(is.null(ffx$Gtc)){ ## contraints if user don't provide them
        mm <- diag(1,nt,nt); mm[lower.tri(mm)] <- 0; #mm[upper.tri(mm)] <- 2
        gesIf[[u]] <- mm
      }else{gesIf[[u]] <- rep(list(ffx$Gtc),length(ffx$Z))}
    }else{ ## if is a normal term
      
      if(purefixedtermss[u] == "(Intercept)"){
        xpp <- as.matrix(baseX[,"(Intercept)"])
        colnames(xpp) <- "(Intercept)"
      }else{
        if(types[purefixedtermss[u]] == "numeric"){ #if is numeric just look for the name of the column
          findlevs <- purefixedtermss[u]
          xpp <- as.matrix(baseX[,findlevs]); colnames(xpp) <- purefixedtermss[u]
        }else{
          findlevs <- paste0( purefixedtermss[u] , unique(as.character(na.omit(data[,purefixedtermss[u]]))) )
          xpp <- as.matrix(baseX[, which(colnames(baseX) %in% findlevs)])
        }
      }
      # print(head(xpp))
      xs[[u]] <- list(xpp)
      gesf[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)) 
      mm <- diag(1,nt,nt); mm[lower.tri(mm)] <- 0; #mm[upper.tri(mm)] <- 2
      gesIf[[u]] <- list(mm)
    }
  }
  # print(str(xs))
  # print(str(gesIf))
  gesIx <- unlist(gesIf,recursive=FALSE)
  X <- unlist(xs,recursive=FALSE)
  ## add intercept
  # intercheck <- which(colnames(baseX) %in% "(Intercept)")
  # if(length(intercheck) >0){
  #   gesIx <- c(list(diag(nt)),gesIx)
  #   px <- as.matrix(baseX[,"(Intercept)"]); colnames(px) <- "(Intercept)"
  #   X <- c(list(px),X)
  # }
  Gx <- gesIx
  # print(str(X))
  # print(str(Gx))
  #################
  ## weights
  # if(missing(weights)){
  #   ws <- rep(1,nrow(yvar))
  # }else{ws <- weights}
  
  if(!missing(weights)){
    col1 <- deparse(substitute(weights))
    coco <- data[[col1]]
    ws<- coco
  }else{ws <- rep(1,nrow(yvar))}
  
  #################
  ## subset data
  if(method == "NR"){
    selected <- FALSE
  }else{selected <- TRUE}
  #################
  ## provide arguments to MNR
  if(!missing(random)){
    Z <- lapply(Z,function(x){as(x, Class = "sparseMatrix")})
  }else{Z <- list();K <- list();random=NULL}
  R <- lapply(R,function(x){as(x, Class = "sparseMatrix")})
  if(!is.null(init)){GES <- init}
  if(!is.null(constraints)){GESI <- constraints}
  re_names <- c(re_namel1,re_namel2)
  
  if(return.param){
    good <- provdat$good
    res <- list(yvar, X,Gx,Z,K,R,GES,GESI, ws,
                iters, tolpar, tolparinv, 
                selected,getPEV,verbose,re_names,
                good,fixedtermss
                )
  }else{
    res <- .Call("_sommer_MNR",PACKAGE = "sommer",yvar, X,Gx,Z,K,R,GES,GESI, ws,
                   iters, tolpar, tolparinv,
                   selected,getPEV,verbose, FALSE)
    
    # res <- MNR(yvar, X,Gx,Z,K,R,GES,GESI, ws,
    #              iters, tolpar, tolparinv,
    #              selected,getPEV,verbose, FALSE)
    
    nslices <- dim(res$sigma)[3]
    itraits <- colnames(yvar)
    re_names <- gsub("\\(Intercept):","",re_names)
    re_names_onlyrandom <- gsub("\\(Intercept):","",re_namel1)
    dimnames(res$sigma) <- list(itraits,itraits,re_names)
    names(res$U) <- re_names_onlyrandom
    names(res$VarU) <- re_names_onlyrandom
    names(res$PevU) <- re_names_onlyrandom
    res$method <- method
    res$call <- list(fixed=fixed,random=random,rcov=rcov,
                     na.method.Y=na.method.Y,na.method.X=na.method.X)
    if(!missing(random)){
      namelist <- lapply(Z,function(x){colnames(x)})
    }else{namelist <- list()}
    # print(Gx)
    namesbeta <- lapply(as.list(1:length(X)),function(x){
      tt <- colnames(yvar)[as.logical(apply(Gx[[x]],1,sum))]
      ttp <- expand.grid(tt,colnames(X[[x]]))
      return(ttp)
    })
    # print(namesbeta)
    namesbeta <- do.call(rbind,namesbeta); 
    namelist[[length(namelist)+1]] <- namesbeta
    # res$namelist <- namelist
    if(reshape.output){
      res <- reshape_mmer(res,namelist) 
    }
    res$constraints <- GESI
    res$constraintsF <- Gx
    res$data <- data#dataor
    class(res)<-c("mmer")
  }

  return(res)
}
