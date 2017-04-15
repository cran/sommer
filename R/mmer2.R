mmer2 <- function(fixed, random, data, G=NULL, W=NULL, method="NR", REML=TRUE, MVM=FALSE, 
                  iters=20, draw=FALSE, init=NULL, family=gaussian, silent=FALSE, constraint=TRUE, 
                  sherman=FALSE, EIGEND=FALSE, forced=NULL, map=NULL, fdr.level=0.05, manh.col=NULL, 
                  min.n=FALSE, gwas.plots=TRUE, n.cores=1, tolpar = 1e-06, tolparinv = 1e-06, 
                  IMP=TRUE, n.PC=0, P3D=TRUE, models="additive", ploidy=2, min.MAF=0.05){
  #rcov <- missing
  #gss=TRUE
  if(missing(data)){
    data <- environment(fixed)
    data2 <- environment(random)
    nodata <-TRUE
    #cat("data argument not provided \n")
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
  
  if(!is.null(G) & method == "EM"){
    cat("With var-cov structures (G) present you may want to try the AI or NR algorithm.\n\n")
  }
  ###########################
  ########### useful functions
  at <- function(x, levs){ # how you handle x
    dd <- model.matrix(~x - 1,data.frame(x))
    colnames(dd) <- substring(colnames(dd),2)
    dd <- dd[,levs]
    return(dd)
  }
  diagc <- grep("diag\\(", random)
  if(length(diagc)){ # if user fits a diagonal model is the same than at
    random <- as.formula(paste(gsub("diag","at",random),collapse = ""))
  }
  g <- function(x){x}
  and <- function(x){x}
  
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  
  #   uscheck <- grep("us\\(", random)
  #   if(length(uscheck)>0){
  #     stop("The us(.) function is not available in sommer yet. \n",call. = FALSE)
  #   }
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
  #print(str(data2))
  mf <- try(model.frame(fixed, data = data2, na.action = na.pass), silent = TRUE)
  mfna <- try(model.frame(fixed, data = data, na.action = na.pass), silent = TRUE)
  if (class(mf) == "try-error") {
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mf <- eval(mf, parent.frame())
  mfna <- eval(mfna, parent.frame())
  
  #which(!duplicated(t(mfna)))
  # response Y
  yvar <- model.response(mfna)
  
  #yvar <- gsub(" ", "", as.character(fixed[2]))
  ### Xb in 'formula'
  #print(str(mf))
  #print(fixed)
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
    #zvar.names <- colnames(V)
    #zvar.names <- zvar.names[which(!duplicated(zvar.names))]
    #print(zvar.names)
    zvar <- V #names(V)
    #     for(i in 1:dim(zvar)[2]){
    #       zvar[,i] <- as.factor(zvar[,i])
    #     }
    #zvar <- apply(zvar,2,as.factor)
    #zvar.names <- names(V)
    #zvar <- gsub(" ", "", strsplit(as.character(random[2]), split = "[+]")[[1]])
    #varsss <- c(xvar,zvar)
    
    ### filters to apply at, us, diag, g
    #     atc <- grep("at",zvar.names)
    #     if(length(atc)>0){
    #       apply(data.frame(zvar.names[atc]),1,function(x){gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]])})
    #       
    #     }
    
    Z <- list()
    counter <- 0
    for(i in 1:length(zvar.names)){
      ## incidence matrix
      vara <- zvar.names[i]
      
      # data.frame(factor(V[,vara],levels=V[,vara],ordered=T))
      zi <- model.matrix(as.formula(paste("~",vara,"-1")),zvar)
      
      ### check for overlay matrices
      andc <- grep("and\\(",vara)
      if(length(andc)>0){
        
        if(!is.factor(zvar[,vara])){stop(paste("Random effect",vara,"needs to be a factor if used as random. \nPlease convert using as.factor() function and fit again.\n"), call.=FALSE)}
        
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
          j3 <- gsub(".*:","",vara) # part removed
          j2 <- gsub(" ","", j2)
          k1 <- gsub(":.*","",vara) # to remove in next step
          #regexpr("\\((.*)\\)", j2[1])
          orx <- gsub(k1,"",j2, fixed = TRUE) # order of columns by location
          where <- as.matrix(apply(data.frame(unique(orx)),1,function(x,y){which(y==x)},y=orx)) # each column says indeces for each level of at()
          colnames(where) <- unique(j2)
          for(u in 1:dim(where)[2]){
            counter <- counter+1
            zix <- zi[,where[,u]]
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
                inte <- intersect(uuuz,uuuk)
                #print(length(inte)==length(uuuz))
                if(length(inte)==length(uuuz)){ # the names were the same in Z and K
                  ki <- G[[ww]][colnames(zix),colnames(zix)]#[uuuz,uuuz]
                }else{ # no intersection between z and k names
                  
                  cat(paste("\nNames of Z and K for random effect",vara,"are not the same. \nMake sure they are in the correct order."))
                  
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
          constraint=FALSE #there can be negative covariance components
          
          j1 <- colnames(zi)
          j2 <- gsub(":.*","",j1) # remove everything after the : to all names
          j3 <- gsub(".*:","",vara) # part removed
          j2 <- gsub(" ","", j2)
          k1 <- gsub(":.*","",vara) # to remove in next step
          #regexpr("\\((.*)\\)", j2[1])
          orx <- gsub(k1,"",j2, fixed = TRUE) # order of columns by location
          where <- as.matrix(apply(data.frame(unique(orx)),1,function(x,y){which(y==x)},y=orx)) # each column says indeces for each level of at()
          colnames(where) <- unique(j2)
          
          df1 <- expand.grid(colnames(where),colnames(where))
          df1 <- df1[!duplicated(t(apply(df1, 1, sort))),]
          rownames(df1) <- NULL
          
          #           if(is.null(init)){
          #           arevar <- which(df1[,1]==df1[,2])
          #           arecovar <- which(df1[,1]!=df1[,2])
          #           lolo<-sum(length(arecovar)+length(arevar))+1
          #           init <- rep(var(yvar, na.rm = TRUE)/lolo,lolo)
          #           init[arecovar] <- (init[arecovar]^2)/1000
          #           print(init)
          #           # first fit a model with no covariance effects and then use those as initial
          #           }
          
          for(u in 1:dim(df1)[1]){
            counter <- counter+1
            
            zz <- where[,df1[u,1]] # columns to take for location i
            zz2 <- where[,df1[u,2]] # columns to take for location j
            
            zix <- zi[,zz] ## Z
            zixt <- zi[,zz2] ## Z'
            
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
                inte <- intersect(uuuz,uuuk)
                #print(length(inte)==length(uuuz))
                if(length(inte)==length(uuuz)){ # the names were the same in Z and K
                  ki <- G[[ww]][colnames(zix),colnames(zix)]#[uuuz,uuuz]
                }else{ # no intersection between z and k names
                  
                  cat(paste("\nNames of Z and K for random effect",vara,"are not the same. \nMake sure they are in the correct order."))
                  
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
            #elem <- list(K=(zix%*%ki%*%t(zixt)))#, Zt=zixt)
            elem <- list(Z=zix, K=ki, Zt=zixt)
            Z[[counter]] <- elem
            names(Z)[counter] <- paste(paste(as.character(unlist(df1[u,])),collapse = ":"),":",j3,sep="")
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
            inte <- intersect(uuuz,uuuk)
            if(length(inte)==length(uuuz)){ # the names were the same in Z and K
              ki <- G[[ww]][colnames(zi),colnames(zi)]#[uuuz,uuuz]
            }else{ # no intersection between z and k names
              cat(paste("\nNames of Z and K for random effect",vara,"are not the same. \nMake sure they are in the correct order."))
              ki <- G[[ww]] 
            }
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
      
    }
    #print(str(Z))
    #names(Z) <- zvar.names
    
    #     if(missing(rcov)){
    #       R <- NULL
    #     }else{
    #       rvar.names <- gsub(" ", "", strsplit(as.character(rcov[2]), split = "[+]")[[1]])
    #       R <- list()
    #       for(n in 1:length(rvar.names)){ #Ri
    #         req.ris <- strsplit(rvar.names[n],":")[[1]] # required Rij's
    #         # 1) ar1, 2) cs, 3) arma
    #         typ.r <- gsub("\\(.*","",req.ris) # type of correlation matrix
    #         typ.v <- gsub(".*\\(","",req.ris); typ.v <- gsub("\\)","",typ.v) # for which variable
    #         Ri <- list() # to store Rij's
    #         rit <- vector(mode="character")# to store type of correlation matrices
    #         for(o in 1:length(typ.r)){#Rij
    #           if(typ.r[o]=="ar1"){
    #             nr <- length(table(data[,typ.v[o]]))
    #             Ri[[o]] <- AR1.mat(.25,nr)
    #             rit[o] <- "AR1"
    #           }else if(typ.r[o]=="cs"){
    #             nr <- length(table(data[,typ.v[o]]))
    #             Ri[[o]] <- CS.mat(.25,nr)
    #             rit[o] <- "CS"
    #           }else if(typ.r[o]=="arma"){
    #             nr <- length(table(data[,typ.v[o]]))
    #             Ri[[o]] <- AR1.mat(.25,nr)
    #             rit[o] <- "AR1"
    #           }else if(typ.r[o]=="id"){
    #             nr <- length(table(data[,typ.v[o]]))
    #             Ri[[o]] <- diag(nr)
    #             rit[o] <- "ID"
    #           }
    #         }
    #         ## once we filled Ri put Ri in R
    #         R[[n]] <- Ri
    #         R[[n]]$type <- rit
    #         #names(R[[n]])[o+1] <- "type"
    #       }## Ri
    #       
    #     } # end for rcov present or not
    #     print(str(Z))
    #     print(str(X))
    #     print(str(yvar))
    res <- mmer(Y=yvar, X=X, Z=Z, W=W, method=method, REML=REML, iters=iters, draw=draw, init=init, 
                silent=silent, constraint=constraint, sherman=sherman, EIGEND=EIGEND, forced=forced, 
                map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=gwas.plots,n.cores=n.cores, 
                MVM=MVM,tolpar = tolpar, tolparinv = tolparinv, IMP=IMP, ploidy=ploidy, models=models)
    
  }else{###only fixed effects
    res <- glm(yvars~X, family=family)
  }
  #########
  return(res)
}
