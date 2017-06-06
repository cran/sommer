mmer2 <- function(fixed, random, rcov, data, G=NULL, W=NULL, method="NR", REML=TRUE, DI=TRUE, MVM=FALSE, 
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
  
  if(!is.null(G) & method == "EM"){
    cat("With dense var-cov structures (G) present you may want to try the AI or NR algorithm.\n")
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
    
    if(!is.null(G)){
      ranused <- zvar.names[grep("g\\(",zvar.names)] # g() random effects
      if(length(ranused)){ # if g() were used 
        diesel <- setdiff(names(G),apply(data.frame(ranused),1,expi)) # did they forget to add a g() and provided the G var-covar matrix?
        if(length(diesel)>0){
          warning(paste("variance-covariance matrices specified in the G argument for:\n",paste(diesel,collapse = ", "),"were not used. Plase use the g() function to use such.\n"), call. = FALSE, immediate. = TRUE)
        }
      }
    }
    
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
      
    }
    #print(str(Z))
    #names(Z) <- zvar.names
    
    if(missing(rcov)){ #### MISSING R ARGUMENT
      R <- NULL
    }else{ ## USER SPECIFIED RCOV ARGUMENT
      
      rvar.names <- gsub(" ", "", strsplit(as.character(rcov[2]), split = "[+]")[[1]])
      
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
    #print(str(R))
    #print(lapply(R,function(x){which(x==1)}))
    
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
    
    
    res <- mmer(Y=yvar, X=X, Z=Z, R=R, W=W, method=method, REML=REML, DI=DI, iters=iters, draw=draw, init=init, 
                silent=silent, constraint=constraint, sherman=sherman, EIGEND=EIGEND, forced=forced, 
                map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=gwas.plots,n.cores=n.cores, 
                MVM=MVM,tolpar = tolpar, tolparinv = tolparinv, IMP=IMP, ploidy=ploidy, models=models, che=FALSE)
    
  }else{###only fixed effects
    res <- glm(yvars~X, family=family)
  }
  #########
  #print(str(Z))
  return(res)
}
