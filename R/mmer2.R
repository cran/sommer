mmer2 <- function(fixed, random, rcov, data, weights, G=NULL, 
                  grouping=NULL, method="NR", init=NULL, iters=20,
                  tolpar = 1e-06, tolparinv = 1e-06, draw=FALSE, 
                  silent=FALSE, constraint=TRUE, forced=NULL, 
                  complete=TRUE, restrained=NULL,na.method.X="exclude",
                  na.method.Y="exclude", REML=TRUE, init.equal=TRUE,
                  return.param=FALSE, date.warning=TRUE){
  #rcov <- missing
  #gss=TRUE
  
  if(!is.null(G)){ # adjust the names of G to avoid errors when using G on overlay() for spacing issues 
    names(G) <- apply(as.data.frame(names(G)),1,function(x){as.character(as.formula(paste("~",x)))[2]})
    #print(names(G))
  }
  
  EIGEND=FALSE
  
  if(missing(data)){
    data <- environment(fixed)
    data2 <- environment(random)
    nodata <-TRUE
    cat("data argument not provided \n")
  }else{nodata=FALSE} 
  
  
  if (!inherits(fixed, "formula")) 
    stop("\nfixed must be a formula")
  if (length(fixed) != 3) 
    stop("\nFixed model formula must be of the form \"resp ~ pred\"")
  if(missing(random)){
    stop("Please use 'lm' for fixed effect models", call. = FALSE)
  } else {
    if (!inherits(random, "formula")) 
      stop("\nrandom must be a formula")
    if (length(random) != 2) 
      stop("\nRandom model formula must be of form \" ~ pred\"")
  }
  
  ###########################
  ## reduce the random formula
  
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  
  spl2D <-  function(x.coord,y.coord,at,at.levels, type="PSANOVA", nseg = c(10,10), pord = c(2,2), degree = c(3,3), nest.div = c(1,1) ) {
    
    interpret.covarrubias.formula <-
      function(formula) {
        env <- environment(formula) 
        if(inherits(formula, "character"))          
          formula <- as.formula(formula)
        tf <- terms.formula(formula, specials = c("SAP", "PSANOVA"))
        terms <- attr(tf, "term.labels")
        nt <- length(terms)
        if(nt != 1)
          stop("Error in the specification of the spatial effect: only a sigle bidimensional function is allowed")
        
        res <- eval(parse(text = terms[1]), envir = env)
        res
      }
    
    bbase <-
      function(X., XL., XR., NDX., BDEG.) {
        # Function for B-spline basis
        dx <- (XR. - XL.)/NDX.
        knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
        P <- outer(X., knots, tpower, BDEG.)
        n <- dim(P)[2]
        D <- diff(diag(n), diff = BDEG. + 1) / (gamma(BDEG. + 1) * dx ^ BDEG.)
        B <- (-1) ^ (BDEG. + 1) * P %*% t(D)
        res <- list(B = B, knots = knots)
        res 
      }
    
    tpower <-
      function(x, t, p) {
        # Function for truncated p-th power function
        return((x - t) ^ p * (x > t))
      }
    
    Rten2 <-
      function(X1,X2) {
        one.1 <- matrix(1,1,ncol(X1))
        one.2 <- matrix(1,1,ncol(X2))
        kronecker(X1,one.2)*kronecker(one.1,X2)
      }
    
    MM.basis <-
      function (x, xl, xr, ndx, bdeg, pord, decom = 1) {
        Bb = bbase(x,xl,xr,ndx,bdeg)
        knots <- Bb$knots
        B = Bb$B
        m = ncol(B)
        n = nrow(B)
        D = diff(diag(m), differences=pord)
        P.svd = svd(crossprod(D))
        U.Z = (P.svd$u)[,1:(m-pord)] # eigenvectors
        d = (P.svd$d)[1:(m-pord)]  # eigenvalues
        Z = B%*%U.Z
        U.X = NULL
        if(decom == 1) {
          U.X = ((P.svd$u)[,-(1:(m-pord))])
          X = B%*%U.X
        } else if (decom == 2){
          X = NULL
          for(i in 0:(pord-1)){
            X = cbind(X,x^i)
          }
        } else if(decom == 3) {
          U.X = NULL
          for(i in 0:(pord-1)){
            U.X = cbind(U.X,knots[-c((1:pord),(length(knots)- pord + 1):length(knots))]^i)
          }
          X = B%*%U.X
        } else if(decom == 4) { # Wood's 2013
          X = B%*%((P.svd$u)[,-(1:(m-pord))])
          id.v <- rep(1, nrow(X))
          D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
          Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
          X <- X%*%Xf
          U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
        }
        list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
      }
    
    ####################
    ### if we want to use at and at.levels
    if(!missing(at)){
      col1 <- deparse(substitute(at))
      dat <- data.frame(x.coord, y.coord, at); colnames(dat) <- c("x.coord","y.coord",col1)
      by <- col1
      if(!missing(at.levels)){
        by.levels=at.levels
      }else{by.levels=NULL}
    }else{
      by=NULL
      by.levels=NULL
      dat <- data.frame(x.coord,y.coord); colnames(dat) <- c("x.coord","y.coord")
    }
    #######################
    
    x.coord <- "x.coord"
    y.coord <- "y.coord"
    
    if(is.null(by)){
      dat$FIELDINST <- "FIELD1"
      by="FIELDINST"
      dat[,by] <- as.factor(dat[,by])
      data0 <- split(dat, dat[,by])
      if(!is.null(by.levels)){
        keep <- names(data0)[which(names(data0) %in% by.levels)]
        
        if(length(keep)==0){stop("The by.levels provided were not found in your dataset.",call. = FALSE)}
        
        data0 <- data0[[keep]]
        if(length(keep)==1){data0 <- list(data0); names(data0) <- keep}
      }
    }else{
      check <- which(colnames(dat)==by)
      if(length(check)==0){stop("by argument not found in the dat provided", call. = FALSE)}else{
        
        missby <- which(is.na(dat[,by]))
        if(length(missby)>0){stop("We will split using the by argument and you have missing values in this column.\nPlease correct.", call. = FALSE)}
        
        dat[,by] <- as.factor(dat[,by])
        data0 <- split(dat, dat[,by])
        
        if(!is.null(by.levels)){
          keep <- names(data0)[which(names(data0) %in% by.levels)]
          
          if(length(keep)==0){stop("The by.levels provided were not found in your dataset.",call. = FALSE)}
          
          data0 <- data0[[keep]]
          if(length(keep)==1){data0 <- list(data0); names(data0) <- keep}
          #print(str(data0))
        }
        
      }
    }
    
    nasx <- which(is.na(dat[,x.coord]))
    nasy <- which(is.na(dat[,y.coord]))
    if(length(nasx) > 0 | length(nasy) >0){
      stop("x.coord and y.coord columns cannot have NA's", call. = FALSE)
    }
    #res <- interpret.covarrubias.formula(formula)
    
    ####
    #### now apply the same to all environments
    multires <- lapply(data0, function(data){
      
      
      x1 <- data[ ,x.coord]
      x2 <- data[ ,y.coord]
      
      #type = type
      
      MM1 = MM.basis(x1, min(x1), max(x1), nseg[1], degree[1], pord[1], 4)
      MM2 = MM.basis(x2, min(x2), max(x2), nseg[2], degree[2], pord[2], 4)
      
      X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
      X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B
      
      c1 = ncol(B1); c2 = ncol(B2)
      
      # Nested bases
      if(nest.div[1] == 1) {
        MM1n <- MM1
        Z1n <- Z1
        c1n <- c1
        d1n <- d1	
      } else {
        MM1n = MM.basis(x1, min(x1), max(x1), nseg[1]/nest.div[1], degree[1], pord[1], 4)
        Z1n <- MM1n$Z
        d1n <- MM1n$d
        c1n <-  ncol(MM1n$B)  					
      }
      if(nest.div[2] == 1) {
        MM2n <- MM2
        Z2n <- Z2
        c2n <- c2
        d2n <- d2	
      } else {
        MM2n = MM.basis(x2, min(x2), max(x2), nseg[2]/nest.div[2], degree[2], pord[2], 4)
        Z2n <- MM2n$Z
        d2n <- MM2n$d
        c2n <-  ncol(MM2n$B)  					
      }
      
      x.fixed <- y.fixed <- ""
      for(i in 0:(pord[1]-1)){
        if(i == 1) 
          x.fixed <- c(x.fixed, x.coord)
        else if( i > 1)
          x.fixed <- c(x.fixed, paste(x.coord, "^", i, sep = ""))
      }
      for(i in 0:(pord[2]-1)){
        if(i == 1) 
          y.fixed <- c(y.fixed, y.coord)
        else if( i > 1)
          y.fixed <- c(y.fixed, paste(y.coord, "^", i, sep = ""))
      }
      xy.fixed <- NULL
      for(i in 1:length(y.fixed)) {
        xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
      }
      xy.fixed <- xy.fixed[xy.fixed != ""]
      names.fixed <- xy.fixed
      
      smooth.comp <- paste("f(", x.coord,",", y.coord,")", sep = "")
      
      if(type == "SAP") {
        names.random <- paste(smooth.comp, c(x.coord, y.coord), sep = "|")				
        X = Rten2(X2, X1)		
        # Delete the intercept
        X <- X[,-1,drop = FALSE]
        Z = cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2n, Z1n))
        
        dim.random <- c((c1 -pord[1])*pord[2] , (c2 - pord[2])*pord[1], (c1n - pord[1])*(c2n - pord[2]))		
        dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
        names(dim$fixed) <- names.fixed
        names(dim$random) <- paste(smooth.comp, "Global")
        
        # Variance/Covariance components
        g1u <- rep(1, pord[2])%x%d1
        g2u <- d2%x%rep(1, pord[1])
        g1b <- rep(1, c2n - pord[2])%x%d1n
        g2b <- d2n%x%rep(1, c1n - pord[1])
        
        g <- list()	
        g[[1]] <- c(g1u, rep(0, dim.random[2]), g1b)
        g[[2]] <- c(rep(0, dim.random[1]), g2u, g2b)
        
        names(g) <- names.random
        
      } else {		
        one1. <- X1[,1, drop = FALSE]
        one2. <- X2[,1, drop = FALSE]
        
        x1. <- X1[,-1, drop = FALSE]
        x2. <- X2[,-1, drop = FALSE]
        
        # Fixed and random matrices
        X <- Rten2(X2, X1)
        # Delete the intercept
        X <- X[,-1,drop = FALSE]
        Z <- cbind(Rten2(one2., Z1), Rten2(Z2, one1.), Rten2(x2., Z1), Rten2(Z2, x1.), Rten2(Z2n, Z1n))
        
        dim.random <- c((c1-pord[1]), (c2-pord[2]), (c1-pord[1])*(pord[2]-1), (c2-pord[2])*(pord[1]-1), (c1n-pord[2])*(c2n-pord[2]))
        
        # Variance/Covariance components		
        g1u <- d1
        g2u <- d2
        
        g1v <- rep(1, pord[2] - 1)%x%d1
        g2v <- d2%x%rep(1,pord[1] - 1)
        
        g1b <- rep(1, c2n - pord[2])%x%d1n
        g2b <- d2n%x%rep(1, c1n - pord[1])
        
        g <- list()
        
        if(type == "SAP.ANOVA") {
          g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
          g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
          g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, dim.random[4]), g1b)
          g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, g2b)
          
          names.random <- c(paste("f(", x.coord,")", sep = ""), paste("f(", y.coord,")", sep = ""), paste(smooth.comp, c(x.coord, y.coord), sep = "|"))			
          dim <- list(fixed = rep(1, ncol(X)), random = c(dim.random[1:2], sum(dim.random[-(1:2)])))		
          names(dim$fixed) <- names.fixed
          names(dim$random) <- c(names.random[1:2], paste(smooth.comp, "Global"))
          names(g) <- names.random
        } else {
          g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
          g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
          g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, sum(dim.random[4:5])))
          g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, rep(0, dim.random[5]))
          g[[5]] <- c(rep(0, sum(dim.random[1:4])), g1b + g2b)
          
          names.random <- c(paste("f(", x.coord,")", sep = ""), paste("f(", y.coord,")", sep = ""),
                            paste("f(", x.coord,"):", y.coord, sep = ""),
                            paste(x.coord,":f(", y.coord,")", sep = ""),
                            paste("f(", x.coord,"):f(", y.coord,")", sep = ""))
          
          dim <- list(fixed = rep(1, ncol(X)), random = dim.random)		
          names(dim$fixed) <- names.fixed
          names(dim$random) <- names.random
          names(g) <- names.random
        }		
      }
      colnames(X) <- names.fixed
      colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
      
      attr(dim$fixed, "random") <- attr(dim$fixed, "sparse") <- rep(FALSE, length(dim$fixed))
      attr(dim$fixed, "spatial") <- rep(TRUE, length(dim$fixed))
      
      attr(dim$random, "random") <- attr(dim$random, "spatial") <- rep(TRUE, length(dim$random)) 
      attr(dim$random, "sparse") <- rep(FALSE, length(dim$random))
      
      terms <- list()
      terms$MM <- list(MM1 = MM1, MM2 = MM2)
      terms$MMn <- list(MM1 = MM1n, MM2 = MM2n)
      #terms$terms.formula <- res
      
      # attr(terms, "term") <- smooth.comp
      
      # Initialize variance components
      init.var <- rep(1, length(g))
      
      res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var)	
      M <- cbind(res$X,res$Z)
      
      return(M)
    })
    
    nrowss <- (unlist(lapply(multires,nrow)))
    nranges <- (unlist(lapply(multires,ncol)))
    
    names(multires) <- gsub(" ",".",names(multires))
    names(multires) <- gsub("#",".",names(multires))
    names(multires) <- gsub("/",".",names(multires))
    names(multires) <- gsub("%",".",names(multires))
    names(multires) <- gsub("\\(",".",names(multires))
    names(multires) <- gsub(")",".",names(multires))
    
    
    st <- 1 # for number of rows
    st2 <- 1 # for number of column
    end2 <- numeric() # for end of number of column
    dataflist <- list()
    #glist <- list()
    for(u in 1:length(multires)){
      prov <- multires[[u]]
      mu <- as.data.frame(matrix(0,nrow = sum(nrowss), ncol = ncol(prov)))
      colnames(mu) <- paste(names(multires)[u],colnames(prov), sep="_")
      
      nam <- paste("at",names(multires)[u],"2Dspl", sep="_")
      end <- as.numeric(unlist(st+(nrowss[u]-1)))
      
      mu[st:end,] <- prov
      dataflist[[nam]] <- as(as.matrix(mu), Class="sparseMatrix")
      st <- end+1
      ## for keeping track of the inits
      # end2 <- as.numeric(unlist(st2+(ncol(prov)-1)))
      # glist[[nam]] <- st2:end2
      # st2 <- end2+1
    }
    
    ## now build the last dataframe and adjust the glist
    # newdatspl <- as.data.frame(do.call(cbind,dataflist))
    # nn <- ncol(dat) # to add to the glist
    # glist <- lapply(glist, function(x){x+nn})
    # newdat <- data.frame(dat,newdatspl)
    
    ## now make the formula
    
    #funny <- paste(paste("grp(",names(dataflist),")",sep=""), collapse=" + ")
    
    ## important
    # newdat: is the neew data frame with original data and splines per location matrices
    # glist: is the argument to provide in group in asreml to indicate where each grouping starts and ends
    # funny: formula to add to your random formula
    dataflist <- lapply(dataflist,as.matrix)
    fin <-dataflist#list(newdat=dataflist, funny=funny) # 
    return(fin)
  }
  
  # see if any of the random terms has eig or group
  #rtermss <-strsplit(as.character(random[2]), split = "[+]")[[1]] #)  gsub(" ", "", 
  
  ################
  #### NA.METHOD.X
  if(na.method.X == "exclude"){
    fufu <- strsplit(as.character(fixed[3]), split = "[+]")[[1]]
    fufu <- apply(data.frame(fufu),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    fufu <- unique(unlist(apply(data.frame(fufu),1,function(x){strsplit(x,":")[[1]]})))
    fufu <- setdiff(fufu,"1")
    if(length(fufu) >0){
      missing1 <- sort(unique(as.vector(unlist(apply(data.frame(data[,fufu]),2,function(x){as.vector(which(is.na(x)))})))))
      if(length(missing1) > 0){
        data <- droplevels(data[-missing1,])
      }
    }
  }else if(na.method.X == "include"){
    fufu <- strsplit(as.character(fixed[3]), split = "[+]")[[1]]
    fufu <- apply(data.frame(fufu),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    fufu <- setdiff(fufu,"1")
    if(length(fufu) >0){
      for(u in fufu){
        v <- which(is.na(data[,u]))
        if(length(v) > 0){ # if there's missing data impute
          #print("imputing")
          data[,u] <- imputev(data[,u])
        }
      }
    }
  }else{
    stop("na.method.X not recognized for sommer", call. = FALSE)
  }
  ##########
  ## NA.METHOD.Y
  if(na.method.Y == "exclude"){
    IMP=FALSE
  }else if(na.method.Y == "include"){
    IMP=TRUE
  }else{
    stop("na.method.Y not recognized for sommer", call. = FALSE)
  }
  
  ######
  yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]]
  rtermss <- apply(data.frame(yuyu),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  
  eigcheck <- grep("eig\\(",rtermss)# apply(data.frame(rtermss),2,function(x){grep("eig",x)})
  if(length(eigcheck)>0){
    EIGEND=TRUE
    if(length(eigcheck)>1){stop("Eigen decomposition is only possible for one random effect.\n", call. = FALSE)}
    f1 <- strsplit(rtermss[eigcheck],":")[[1]]
    f00 <- grep("trait",f1) # position of structure for trait
    f0 <- f1[f00] # structure for eigend
    f2 <- f1[setdiff(1:length(f1),f00)]
    f3 <- paste("g(",expi2(f2),")",sep="")
    if(length(f0)>0){f3 <- paste(f0,f3,sep=":")}
    rtermss <- c(f3,rtermss[-eigcheck])
    random <- as.formula(paste("~",paste(rtermss, collapse = " + "))) # new random
    random <- as.formula(paste(as.character(random),collapse=""))
    #print(random)
  }
  grpcheck <- grep("grp\\(",rtermss)# apply(data.frame(rtermss),2,function(x){grep("eig",x)})
  if(length(grpcheck)>0){
    #
    ZKgrouping.str <- character()
    ZKgrouping <- list()
    cous <- 0
    cous2 <- numeric()
    counter <- 0
    for(h in rev(grpcheck)){ # h <- grpcheck[2]
      counter <- counter+1
      f1 <- strsplit(rtermss[h],":")[[1]]
      #f1 <- apply(data.frame(rtermss[grpcheck]),1,function(x){strsplit(x,":")[[1]]})
      f00 <- grep("trait",f1) # position of structure for trait
      if(length(f00)>0){
        ZKgrouping.str[counter] <- f1[f00]
      }else{
        ZKgrouping.str[counter] <- rep("diag(trait)",length(f1)) 
      }
      f0 <- f1[f00] # structure
      f2 <- f1[setdiff(1:length(f1),f00)] # actual random effect
      #if(length(f0)>0){f2 <- paste(f0,f2,sep=":")}
      random <- as.formula(paste("~",paste(rtermss[-h], collapse = " + "))) # new random
      random <- as.formula(paste(as.character(random),collapse=""))
      namere <- expi2(f2) # name of the random effect to look in to grouping argument
      
      rtermss <- rtermss[-h]
      
      #ZKgrouping.str <- character()
      #for(u in 1:length(namere)){ # u <- 1
      cous <- cous+1
      cous2[counter] <- cous
      grpu <- which(names(grouping)==namere)
      Zgrpu <- grouping[[grpu]]
      if(is.null(Zgrpu)){stop("Random effect specified with the grp() function not specified in the grouping argument.\n",call. = FALSE)}
      Gu <- which(names(G)==namere)
      if(length(Gu) > 0){
        Kgrp <- G[[Gu]]
        cat("Ignore warning messages about variance-covariance matrices specified in the \nG argument not used for your grouping effects. They will be used.\n")
      }else{ Kgrp <- diag(ncol(Zgrpu))}
      
      ZKgrouping[[namere]] <- list(Z=Zgrpu, K=Kgrp)
      #names(ZKgrouping)[counter] <- namere
      #}
      
    }
    counter.grp <- counter
    
  }
  if(length(grpcheck)>0){
    cous2.grp <- cous2
  }
  
  ## we will keep 'counter.grp' which is the last grouping factor (length(grp terms))
  ## also keep 'cous2.grp' which it says which terms are the grouping (i.e. 1,2)
  
  splcheck <- grep("spl2D\\(",rtermss)# apply(data.frame(rtermss),2,function(x){grep("eig",x)})
  if(length(splcheck)>0){
    #
    ZKspl2d.str <- character()
    ZKspl2d <- list()
    cous <- 0
    cous2 <- numeric()
    counter <- 0
    for(h in rev(splcheck)){ # h <- grpcheck[2]
      counter <- counter+1
      f1 <- strsplit(rtermss[h],":")[[1]]
      #f1 <- apply(data.frame(rtermss[grpcheck]),1,function(x){strsplit(x,":")[[1]]})
      f00 <- grep("trait",f1) # position of structure for trait
      f0 <- f1[f00] # structure
      f2 <- f1[setdiff(1:length(f1),f00)] # actual random effect
      
      random <- as.formula(paste("~",paste(rtermss[-h], collapse = " + "))) # new random
      random <- as.formula(paste(as.character(random),collapse=""))
      
      rtermss <- rtermss[-h]
      
      provspl <- eval(parse(text = f2),envir = data)
      provspl <- lapply(provspl, function(x){zxk <- list(Z=x, K=diag(ncol(x))); return(zxk)})
      ZKspl2d <- c(ZKspl2d,provspl)
      
      # if(length(f00)>0){
      #   ZKgrouping.str[counter] <- f1[f00]
      # }else{
      #   ZKgrouping.str[counter] <- rep("diag(trait)",length(f1)) 
      # }
      if(length(f00)>0){
        ZKspl2d.str <- c(ZKspl2d.str,rep(f1[f00],length(provspl)))
      }else{
        ZKspl2d.str <- c(ZKspl2d.str,rep("diag(trait)",length(provspl)) )
      }
      provspl <- NULL
    }
    counter.spl2d <- length(ZKspl2d.str)
    
  }
  
  if(length(splcheck)>0){
    cous2.spl2d <- 1:length(ZKspl2d.str)#cous2
  }
  
  ## for us structures
  # apply(data.frame(rtermss),2,function(x){grep("eig",x)})
  # if(length(usscheck)>0 & constraint){
  #   cat("Setting 'constraint' argument to FALSE for unstructured model. \nPlease be carefult with results",call. = FALSE)
  #   constraint=FALSE
  # }
  
  if(missing(rcov)){rcovi <- NULL}else{rcovi<-rcov}
  ttt <- vctable.help(random = random, rcov = rcovi) # extract trait structure
  random <- ttt$random # update random term
  random <- as.formula(paste(as.character(random),collapse=""))
  
  rcov <- ttt$rcov # update rcov term
  
  
  # change the trait structure for residuals from 'at' to 'diag'
  # no valid for me to make 'at' and assume that traits share the same residual value
  ttt$rcov.trt.str <- gsub("at","diag",ttt$rcov.trt.str)
  
  usscheck <- grep("us\\(",strsplit(as.character(random[2]), split = "[+]")[[1]])
  ###########################
  ########### useful functions
  at <- function(x, levs){ 
    dd <- model.matrix(~x - 1,data.frame(x))
    colnames(dd) <- substring(colnames(dd),2)
    dd <- dd[,levs]
  }
  
  diagc <- grep("diag\\(", random)
  if(length(diagc)){ # if user fits a diagonal model is the same than at
    random <- as.formula(paste(gsub("diag","at",random),collapse = ""))
    random <- as.formula(paste(as.character(random),collapse=""))
  }
  
  if(!missing(rcov)){
    diagr <- grep("diag\\(", rcov)
    if(length(diagr)){ # if user fits a diagonal model is the same than at
      rcov <- as.formula(paste(gsub("diag","at",rcov),collapse = ""))
    }
  }
  
  us <- function(x){x}
  g <- function(x){x}
  and <- function(x){x}
  eig <- function(x){x}
  perc.na <- function(x){length(which(is.na(x)))/length(x)}
  
  ###########################
  ############## impute data
  ginvcheck <- grep("ginv\\(",random)
  if(length(ginvcheck)>0){
    #?mmer2
    stop("The homolog to the 'ginv(.)' function in sommer is used as 'g(.)' and the \ncovariance matrices are specifid in the 'G' argument (not the inverse). ",call. = FALSE)
  }
  
  data2 <- data
  if (length(dim(data2)[2]) == 0) {
    stop("Please provide the 'data' argument.\n", call. = FALSE)
  }
  for(i in 1:dim(data2)[2]){
    classdata <- unlist(lapply(data2,class))
    x <- data2[,i]
    if(classdata[i]=="numeric"){
      isNA <- which(is.na(x))
      if(length(isNA) >0 & perc.na(x) < 1){x[isNA] <- mean(x,na.rm=TRUE)}
    }else if(classdata[i]=="factor"){
      isNA <- which(is.na(x))
      if(length(isNA) >0 & perc.na(x) < 1){x[isNA] <- as.factor(names(which(table(x) == max(table(x)))[1]))}
    }else if(classdata[i]=="character"){
      isNA <- which(is.na(x))
      if(length(isNA) >0 & perc.na(x) < 1){x[isNA] <- as.character(names(which(table(x) == max(table(x)))[1]))}
    }
    data2[,i] <- x
  }
  ######################
  #print(str(data2))
  mf <- try(model.frame(fixed, data = data2, na.action = na.pass), silent = TRUE)
  mfna <- try(model.frame(fixed, data = data, na.action = na.pass), silent = TRUE)
  if (class(mf) == "try-error") {
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mf <- eval(mf, parent.frame())
  mfna <- eval(mfna, parent.frame())
  
  yvar <- model.response(mfna)
  #print(str(yvar))
  
  X <- model.matrix(fixed, mf)
  
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
    
    yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    zvar.names <- apply(data.frame(yuyu),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    #zvar.names <- strsplit(as.character(random[2]), split = "[+]")[[1]] #) gsub(" ", "", 
    
    if(!is.null(G)){
      ranused <- zvar.names[grep("g\\(",zvar.names)] # g() random effects
      if(length(ranused)){ # if g() were used 
        diesel <- setdiff(names(G),apply(data.frame(ranused),1,expi)) # did they forget to add a g() and provided the G var-covar matrix?
        if(length(diesel)>0){
          warning(paste("variance-covariance matrices specified in the G argument for:\n",paste(diesel,collapse = ", "),"were not used. Plase use the g() function to use such.\n"), call. = FALSE, immediate. = TRUE)
        }
      }
    }
    
    zvar <- V #names(V)
    
    Z <- list()
    counter <- 0
    counterl <- numeric() # to store the names of the random effects for each random effect specified by
    #for example it specifies where each random effects starts, 
    # if at(LOC):x and loc has 2 levels then counterl[i] <- 2
    # so ends up being something like 2 4 6
    
    # users, specially for interaction where a single term (i.e. at(location):gca) can produce several
    for(i in 1:length(zvar.names)){
      ## incidence matrix
      vara <- zvar.names[i]
      
      # data.frame(factor(V[,vara],levels=V[,vara],ordered=T))
      zi <- model.matrix(as.formula(paste("~",vara,"-1")),zvar)
      
      ### check for overlay matrices
      andc <- grep("and\\(",vara)
      if(length(andc)>0){
        
        #as.formula(paste("~",vara))
        #if(!is.factor(zvar[,vara])){stop(paste("Random effect",vara,"needs to be a factor if used as random. \nPlease convert using as.factor() function and fit again.\n"), call.=FALSE)}
        
        if(i==1){
          stop("To use the overlay a random effect must preceed the and(.) function.", call. = FALSE)
        }
        toverlay <- c(zvar.names[i-1],zvar.names[i])
        xx <- zvar[,toverlay]
        ### check if previous random effect had a g structure to adjust names
        if(length(grep("g\\(",zvar.names[i-1]))>0){ # adjust both names
          colnames(xx)[1:2] <- c(expi(colnames(xx)[1]),expi(colnames(xx)[2]))
        }else{#only adjust the current random effect name
          colnames(xx)[2] <- expi(colnames(xx)[2])
        }
        ##
        zi <- overlay(xx)
        counter <- counter - 1 # move back the counter to replace the model matrix
      }
      ### part to store the matrices
      #### zix <- zi[,df1[u,1]]
      atc <- grep("at\\(",vara)
      usc <- grep("us\\(",vara)
      if(length(atc)>0 | length(usc)>0){
        if(length(atc)>0){ ##$$ IF USER HAS at()$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          if(length(diagc)>0){ #if user wanted diagonal we change names of at to diag
            colnames(zi) <- gsub("at\\(","diag\\(",colnames(zi))
          }
          j1 <- colnames(zi)
          j2 <- gsub(":.*","",j1) # remove everything after the : to all names
          
          #?gregexpr
          st1 <- gregexpr(":",vara)[[1]][1]+1
          j3 <- substring(vara,st1, nchar(vara))
          
          #j3 <- gsub(".*:","",vara) # part removed, should be Row
          #j2 <- gsub(" ","", j2)
          k1 <- gsub(":.*","",vara) # to remove in next step
          #regexpr("\\((.*)\\)", j2[1])
          orx <- gsub(k1,"",j2, fixed = TRUE) # order of columns by location
          where <- matrix(apply(data.frame(unique(orx)),
                                1,function(x,y){which(y==x)},y=orx),
                          ncol=length(unique(orx))) # each column says indeces for each level of at()
          colnames(where) <- unique(j2)
          for(u in 1:dim(where)[2]){ # u=1
            counter <- counter+1
            zix <- as.matrix(zi[,where[,u]])
            ##
            ## var-cov matrix
            gcheck <- grep("g\\(",vara)
            if(length(gcheck)>0){
              if(is.null(G)){stop("You have specified a 'g(.)' structure but G parameter is null\n", call. = FALSE)}
              ww <- which(names(G) %in% expi(vara))
              if(length(ww) > 0){# K was provided
                ## just if there's a K matrix we make sure to be using the real names and no the model.matrix ones
                zas <- gsub(paste("g\\(",expi(vara)[2],")",sep=""),"",colnames(zix))# we remove the g(.) text
                ox <- gsub(",",", ", colnames(where)[u])
                ox <- gsub("\\(","\\\\(",ox)
                colnames(zix) <- gsub(paste(ox,":",sep=""),"",zas)
                #gsub(paste("at\\",substring(zas2,3),":",sep=""),"",zas)
                #colnames(zix) <- gsub(paste("at\\",substring(colnames(where)[u],3),":",sep=""),"",zas) # we remove the at(.) text
                #########
                uuuz <- colnames(zix)#levels(as.factor(colnames(zi))) # order of Z
                uuuk <- attr(G[[ww]],"dimnames")[[1]] # order of K
                if(is.null(uuuk)){uuuk <- dimnames(G[[ww]])[[1]]}
                inte <- intersect(uuuz,uuuk)
                #print(length(inte)==length(uuuz))
                if(length(inte)==length(uuuz)){ # the names were the same in Z and K
                  ki <- G[[ww]][colnames(zix),colnames(zix)]#[uuuz,uuuz]
                }else{ # no intersection between z and k names
                  
                  cat(paste("\nLevels for the random effect",vara," do not coincide with the levels in it variance covariance matrix (G argument).\n"))
                  
                  ki <- G[[ww]] 
                }
                
              }else{ # g was indicated but was not provided, we create a diagonal
                if(u==1){
                  cat(paste("'g(.)' structure indicated for",vara,"but not found in the G parameter.\n"))
                }
                ki <- diag(dim(zix)[2])
              }
            }else{ # was not provided, we create a diagonal
              ki <- diag(dim(zix)[2])
            }
            ###
            elem <- list(Z=zix, K=ki)
            Z[[counter]] <- elem
            
            names(Z)[counter] <- paste(colnames(where)[u],":",j3,sep="")
          }
        }
        #### end of if atc()
        if(length(usc)>0){ ##$$ IF USER HAS us()$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          #constraint=FALSE #there can be negative covariance components
          
          j1 <- colnames(zi)
          j2 <- gsub(":.*","",j1) # remove everything after the : to all names
          j3 <- gsub(".*:","",vara) # part removed
          #j2 <- gsub(" ","", j2)
          k1 <- gsub(":.*","",vara) # to remove in next step
          #regexpr("\\((.*)\\)", j2[1])
          orx <- gsub(k1,"",j2, fixed = TRUE) # order of columns by location
          where <- as.matrix(apply(data.frame(unique(orx)),1,function(x,y){which(y==x)},y=orx)) # each column says indeces for each level of at()
          colnames(where) <- unique(j2)
          
          df1 <- expand.grid(colnames(where),colnames(where))
          df1 <- df1[!duplicated(t(apply(df1, 1, sort))),]
          rownames(df1) <- NULL
          
          for(u in 1:dim(df1)[1]){
            counter <- counter+1
            
            zz <- where[,df1[u,1]] # columns to take for location i
            zz2 <- where[,df1[u,2]] # columns to take for location j
            
            zix <- as.matrix(zi[,zz]) ## Z
            zixt <- as.matrix(zi[,zz2]) ## Z'
            
            ##
            ## var-cov matrix
            gcheck <- grep("g\\(",vara)
            if(length(gcheck)>0){
              if(is.null(G)){stop("You have specified a 'g(.)' structure but G parameter is null\n", call. = FALSE)}
              ww <- which(names(G) %in% expi(vara))
              if(length(ww) > 0){# K was provided
                ## just if there's a K matrix we make sure to be using the real names and no the model.matrix ones
                zas <- gsub(paste("g\\(",expi(vara)[2],")",sep=""),"",colnames(zix))# we remove the g(.) text
                ox <- gsub(",",", ", df1[u,1])
                ox <- gsub("\\(","\\\\(",ox)
                colnames(zix) <- gsub(paste(ox,":",sep=""),"",zas)
                
                zas2 <- gsub(paste("g\\(",expi(vara)[2],")",sep=""),"",colnames(zixt))# we remove the g(.) text
                ox2 <- gsub(",",", ", df1[u,2])
                ox2 <- gsub("\\(","\\\\(",ox2)
                colnames(zixt) <- gsub(paste(ox2,":",sep=""),"",zas2)
                
                ######### just one is needed since z and z' have the same g
                uuuz <- colnames(zix)#levels(as.factor(colnames(zi))) # order of Z
                uuuk <- attr(G[[ww]],"dimnames")[[1]] # order of K
                if(is.null(uuuk)){uuuk <- dimnames(G[[ww]])[[1]]}
                inte <- intersect(uuuz,uuuk)
                #print(length(inte)==length(uuuz))
                if(length(inte)==length(uuuz)){ # the names were the same in Z and K
                  ki <- G[[ww]][colnames(zix),colnames(zix)]#[uuuz,uuuz]
                }else{ # no intersection between z and k names
                  
                  cat(paste("\nNames of Z and K for random effect",vara,"do not all match. Make sure they are in the correct order.\n"))
                  
                  ki <- G[[ww]] 
                } # end of if g was provided
                
              }else{ # g was indicated but was not provided, we create a diagonal
                if(u==1){
                  cat(paste("'g(.)' structure indicated for",vara,"but not found in the G parameter.\n"))
                }
                ki <- diag(dim(zix)[2])
              }
            }else{ # was not provided, we create a diagonal
              ki <- diag(dim(zix)[2])
            }
            ###
            vc.check <- as.character(unlist(df1[u,]))
            if(vc.check[1] == vc.check[2]){# is variance
              elem <- list(Z=zix, K=ki)
            }else{# is covariance
              ki0 <- ki*0
              zix.cov <- cbind(zix,zixt)
              ki.cov <- rbind(cbind(ki0,ki),cbind(ki,ki0))
              elem <- list(Z=zix.cov, K=ki.cov)
            }
            
            #elem <- list(Z=zix, K=ki)
            Z[[counter]] <- elem
            # Env:Name!ca:fl
            env <- expi2(as.character(unlist(df1[u,])))[1]# 
            env.name <- paste(env,j3,sep=":") # Env
            inter <- gsub(paste("us\\(",env,")",sep=""),"",paste(as.character(unlist(df1[u,])),collapse = ":"))
            names(Z)[counter] <- paste(env.name,paste(env,inter,sep="."),sep="!")
            
            #names(Z)[counter] <- paste(paste(as.character(unlist(df1[u,])),collapse = ":"),":",j3,sep="")
          }
          ##specify inital values if unstructured model
          
          ## end of specifying initial values
        }
        #### end of if usc()
      }else{ ##$$$$$$ if user is not using at() or us() $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        counter <- counter+1
        ## var-cov matrix
        gcheck <- grep("g\\(",vara)
        if(length(gcheck)>0){
          if(is.null(G)){stop("You have specified a 'g(.)' structure but G parameter is null\n", call. = FALSE)}
          ww <- which(names(G) %in% expi(vara))
          #% if and(g(.))
          if(length(andc)>0){ww <- which(names(G) %in% expi2(vara))}
          #$
          if(length(ww) > 0){# K was provided
            ## just if there's a K matrix we make sure to be using the real names and no the model.matrix ones
            colnames(zi) <- gsub(paste("g\\(",expi(vara),")",sep=""),"",colnames(zi))#levels(as.factor(V[,vara]))
            
            #########
            uuuz <- colnames(zi)#levels(as.factor(colnames(zi))) # order of Z
            uuuk <- attr(G[[ww]],"dimnames")[[1]] # order of K
            if(is.null(uuuk)){uuuk <- dimnames(G[[ww]])[[1]]}
            #print(uuuk)
            ### ========== sommer 2.8 =========== ###
            ## no intersection when "Year:" is present
            dotcheck <- grep(":g\\(",vara) # if user has Year:g(id)
            if(length(dotcheck)>0){
              
              ziplist <- list()
              kiplist <- list()
              uuuc <- unique(gsub(":.*","",uuuz)) # number of levesl before :
              for(h in 1:length(uuuc)){ # now reform Z and form K
                hehe <- grep(uuuc[h],uuuz)
                tostore0 <- uuuz[hehe]
                # form Z 
                ziplist[[h]] <- zi[,tostore0]
                # now K
                tostore1 <- gsub(".*:","",tostore0)
                kipprov <- G[[ww]][tostore1,tostore1]
                colnames(kipprov) <- rownames(kipprov) <- tostore0
                kiplist[[h]] <- kipprov
              }
              
              zi <- do.call("cbind",ziplist)
              ki <- do.call("adiag1",kiplist)
              
            }else{
              ##$$
              inte <- intersect(uuuz,uuuk)
              if(length(inte)==length(uuuz)){ # the names were the same in Z and K
                ki <- G[[ww]][colnames(zi),colnames(zi)]#[uuuz,uuuz]
              }else{ # no intersection between z and k names
                cat(paste("\nNames of Z and K for random effect",vara,"are not the same. \nMake sure they are in the correct order."))
                ki <- G[[ww]] 
              }
              ##$$
            }
            ### ================================== ###
            
            
          }else{ # g was indicated but was not provided, we create a diagonal
            cat(paste("'g(.)' structure indicated for",vara,"but not found in the G parameter.\n"))
            ki <- diag(dim(zi)[2])
          }
        }else{ # was not provided, we create a diagonal
          ki <- diag(dim(zi)[2])
        }
        elem <- list(Z=zi, K=ki)
        Z[[counter]] <- elem
        
        names(Z)[counter] <- vara # add the name of the random effect
      }
      counterl[i] <- counter
    }
    
    ##################################
    ####now add the grouping matrices to the model
    ##################################
    if(length(grpcheck)>0){
      Z <- c(Z,ZKgrouping)
      cous2 <- max(counter) + counter.grp # total number of ZK matrices added, counter is only from RE and counter.grp from grouping
      counterl <- c(counterl,max(counter)+cous2.grp) # 
      # originally counterl says the last effect where the same trait structure should be applied, i.e.
      # for at(E):B + at(F):G if 'E' has 3 levels and 'F' 4 levels counterl is c(3,7)
      ttt$ran.trt.str <- c(ttt$ran.trt.str,ZKgrouping.str) 
      # the trait structures should match the length of counterl
    }
    
    if(length(splcheck)>0){
      Z <- c(Z,ZKspl2d)
      if(length(cous2)==0){cous2 <- max(counter)}
      cous2 <- cous2 + counter.spl2d # total number of ZK matrices added
      counterl <- c(counterl,(cous2-counter.spl2d)+cous2.spl2d) # which have a unique trait structure
      ttt$ran.trt.str <- c(ttt$ran.trt.str,ZKspl2d.str)
    }
    
    if(missing(rcov)){ #### MISSING R ARGUMENT
      R <- NULL
    }else{ ## USER SPECIFIED RCOV ARGUMENT
      
      rvar.names <- strsplit(as.character(rcov[2]), split = "[+]")[[1]]#) gsub(" ", "", 
      
      ## control for only units on the meantime
      uni <- gsub(".*:","",rvar.names)
      if(uni != "units"){
        stop("On the meantime the only rcov structures available are:\n 'rcov=~units' or 'rcov=~at(.):units'.",call. = FALSE)
      }
      atrcov <- grep("at\\(",rvar.names) #check
      #if "at" is present:
      if(length(atrcov)>0){
        atroz <- expi2(rvar.names) # variable where we defined the at
        ri <- model.matrix(as.formula(paste("~",atroz,"-1")),data)
        R <- list()
        uuur <- as.character(unique(data[,atroz]))#unique(gsub(":.*","",colnames(ri))) # number of levesl before :
        for(h in 1:length(uuur)){ # now reform Z and form K
          rprov <- matrix(0,dim(data)[1],dim(data)[1])
          year <- uuur[h]#gsub(atroz,"",uuur[h]) # level to look for
          tokeepp <- which(data[,atroz] == year)
          diag(rprov)[tokeepp] <- 1
          # form R
          R[[h]] <- rprov
        }
        #print(uuur)
        names(R) <- paste(uuur,":units",sep="")
        
      }else{#if no at present just units
        R <- NULL
      }
      
    }
    
    ### applying check beffore fitting ===============
    cococo <- names(Z)
    provided <- lapply(Z, names)
    for(s in 1:length(provided)){ #for each random effect =============================
      provided2 <- names(Z[[s]])
      #----the 's' random effect has two matrices----
      dido<-lapply(Z[[s]], dim) # dimensions of Z and K
      condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
      # condition, column size on Z matches with a square matrix K
      if(!condi){
        cat(paste("In the",cococo[s],"random effect you have an incidence matrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"effect levels. \nTherefore the variance-covariance matrix for this random effect in the G argument should be a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2],"matrix. Please check your matrices in G.")
        stop()
      }
    } #for each random effect end =================================================
    
    yvar <- as.matrix(yvar); 
    if(ncol(yvar) == 1){ # to make sure that even for single trait we get the colnames
      dog <- as.character(fixed[[2]]); 
      dog <- setdiff(dog,"cbind")
      colnames(yvar) <- dog
    }
    
    ## do the parameter restrain based on trait structures made at the beggining
    rs <- length(Z) # random effects
    rrs <- length(R) # residual random effects
    if(rrs==0){rrs <- 1}
    mapping <- as.data.frame(vctable(ncol(yvar), rs, rrs))
    #mapping$stru <- NA
    
    strmapping <- data.frame(x=c(counterl[1],diff(counterl),rrs),y=c(ttt$ran.trt.str,ttt$rcov.trt.str))
    ttt.all <- as.vector(unlist(apply(strmapping,1,function(x){rep(x[2],x[1])})))
    #     for(m in 1:length(strmapping)){
    #       mapping$stru[which(mapping$res == m)] <- strmapping[m]
    #     }
    
    #ttt.all <- c(ttt$ran.trt.str, ttt$rcov.trt.str) # structures for each random effect
    torestrain <- numeric()
    for(f in 1:length(ttt.all)){
      structure <- gsub("\\(.*","",ttt.all[f]) # trait structure of the RE
      if(structure == "diag"){
        torestrain <- c(torestrain,which(mapping$t1 != mapping$t2 & mapping$res == f))
      }else if(structure == "us"){
        # if us there's nothing to restrain
      }else{
        stop("Structure for trait not recognized", call. = FALSE)
      }
    }
    if(length(torestrain)==0){torestrain <- NULL}
    
    #print(names(Z))
    if(length(usscheck)>0){ # provide good initial values
      termos <- gsub(".*:","",rtermss[usscheck]) #obtain the re where us is being applied
      ussnames <- apply(data.frame(termos),1,function(x){gsub(paste(":",x,sep=""),"",names(Z))}) # names without the terms
      covars <- apply(data.frame(ussnames),1,function(x){uuu <- strsplit(x,":")[[1]]; return(uuu[1]==uuu[2])})# TRUES are variances, FALSE are covariances
      inicio <- rep(list(var(yvar,na.rm = TRUE)/(length(covars)*.5)),rs+rrs)
      toreduce <- which(!covars) # covariance components need to be reduced
      for(k in 1:length(toreduce)){
        inicio[[toreduce[k]]] <- inicio[[toreduce[k]]]/3
      }
      init <- inicio
    }else{
      termos <- gsub(".*:","",names(Z)) #remove everything before the 1st ':'
      
      rtermss2 <- apply(data.frame(rtermss),1,function(x,yy){
        for(o in 1:length(yy)){
          x <- gsub(yy[o],"",x)
        }
        return(x)
      },yy=unique(termos))
      rtermss2 <-gsub(":","",rtermss2)
      rtermss3 <-gsub("diag\\(","diag\\\\(",rtermss2)
      rtermss4 <-gsub("at\\(","at\\\\(",rtermss2) # form1
      rtermss5 <-gsub("diag\\(","at\\\\(",rtermss2) # form2
      rtermss6 <-gsub("at\\(","diag\\\\(",rtermss2) #form3
      rtermss7 <-c(rtermss3,rtermss4,rtermss5,rtermss6)
      rtermss7 <-trimws(rtermss7)
      
      expi3 <- function(j){gsub("[diag\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
      expi4 <- function(x){gsub("(?<=diag\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
      
      expi3(names(Z))
      
      apply(data.frame(names(Z)),1,function(x,yy){
        for(o in 1:length(yy)){
          x <- gsub(yy[o],"",x)
        }
        return(x)
      },yy=rtermss7)
    }
    #print(torestrain)
    #print(str(Z))
    
    vara2 <- gsub( " *diag\\(.*?\\)) *", "", names(Z))
    vara2 <- gsub( " *at\\(.*?\\)) *", "", vara2)
    vara2 <- gsub( " *us\\(.*?\\)) *", "", vara2)
    vara2 <- gsub("diag\\s*\\([^\\)]+\\)","",vara2)
    vara2 <- gsub("at\\s*\\([^\\)]+\\)","",vara2)
    vara2 <- gsub("us\\s*\\([^\\)]+\\)","",vara2)
    names(Z) <- vara2
    
    if(!missing(weights)){
      
      col1 <- deparse(substitute(weights))
      coco <- data[[col1]]
      W <- Diagonal(n=length(coco),x=coco)
      coco <- NULL; col1 <- NULL
    }else{W<- NULL}
    
    ## force the covariance matrices for overlay
    overc <- grep("overlay\\(",names(Z)) # overlay check
    if(length(overc) > 0){
      for(oo in 1:length(overc)){
        nameoverlay <- names(Z)[overc[oo]] # name of the overlay random effect
        #print(nameoverlay)
        #print(names(G))
        if(!is.null(G)){ # if there's Gs look for overlay Gs
          gover <- grep(gsub("overlay\\(","overlay\\\\(",names(G)), nameoverlay)
        }else{
          gover <- numeric()
        }
        
        ## if user specified the g() for an overlay effect but there was no match
        # if(length(gover)==0){
        #   warning(paste("covariance matrices in G not found for",nameoverlay,". covariance matrices available are:", paste(names(G), collapse = ", ")),call. = FALSE)
        # }
        
        #print(gover)
        if(length(gover)>0){# if there's a covariance structure for the over lay
          Kprov <- G[[gover]]
          if(is.null(colnames(Kprov))){
            stop("Column and row names for the covariance structures are mandatory \nand should match the level names of the random effect.\n", call. = FALSE)
            #colnames(Kprov) <- rownames(Kprov) <- 1:nrow(Kprov)
            }
          #print(Kprov)
          Zprov <- Z[[overc[oo]]]$Z
          correct <- gsub("overlay\\(","overlay\\\\(",nameoverlay)
          
          #print(correct)
          # correct <- gsub("at\\(","at\\\\(",correct)
          # correct <- gsub("diag\\(","diag\\\\(",correct)
          # correct <- gsub("us\\(","us\\\\(",correct)
          colnames(Zprov) <- gsub(paste(".*",correct,sep=""),"",colnames(Zprov))
          #print(colnames(Zprov))
          #print(nameoverlay)
          #print(zprov)
          #colnames(Kprov) <- rownames(Kprov) <- paste(nameoverlay,rownames(Kprov),sep="")
          #print(Kprov)
          #print(Zprov)
          #print(colnames(zprov))
          #Kprov <- Kprov[colnames(zprov),colnames(zprov)]
          xoxo <- overc[oo]
          #print(Kprov)
          Kprov <- Kprov[colnames(Zprov),colnames(Zprov)]
          #print(Kprov)
          Z[[xoxo]]$K <- Kprov #Kprov#G[[gover]]
          Z[[xoxo]]$Z <- Zprov
          #print(str(Z[[xoxo]]))
          #print(overc[oo])
        }
      }
      #print(str(Z))
    }
    #print(str(Z))
    #print(str(Z))
    if(return.param){
      res <- list(Y=yvar, X=X, Z=Z, R=R, W=W, method=method, REML=REML, init=init,
                  iters=iters,tolpar=tolpar,
                  tolparinv = tolparinv,draw=draw,silent=silent, 
                  constraint = constraint,EIGEND = EIGEND,
                  forced=forced,IMP=IMP,complete=complete,restrained=torestrain,
                  init.equal = init.equal)
    }else{
    res <- mmer(Y=yvar, X=X, Z=Z, R=R, W=W, method=method, REML=REML, init=init,
                iters=iters,tolpar=tolpar,
                tolparinv = tolparinv,draw=draw,silent=silent, 
                constraint = constraint,EIGEND = EIGEND,
                forced=forced,IMP=IMP,complete=complete,restrained=torestrain,
                init.equal = init.equal, date.warning=date.warning)
    if(!is.null(res$restrained)){
      names(res$restrained) <- rownames(res$sigma.scaled)[res$restrained] 
    }
    }
    res$call$fixed <- fixed
    res$call$random <- random
    res$call$rcov <- rcov
    res$data <- data
  }else{###only fixed effects
    res <- glm(yvars~X)
  }
  #########
  #print(str(Z))
  return(res)
}
