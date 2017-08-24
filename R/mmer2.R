mmer2 <- function(fixed, random, rcov, data, G=NULL, grouping=NULL, method="NR",
                  init=NULL,iters=20,tolpar = 1e-06, tolparinv = 1e-06, 
                  draw=FALSE, silent=FALSE,
                  constraint=TRUE, EIGEND=FALSE,
                  forced=NULL,IMP=FALSE, complete=TRUE, restrained=NULL){
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

  ###########################
  ## reduce the random formula
  
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  
  # see if any of the random terms has eig or group
  #rtermss <-strsplit(as.character(random[2]), split = "[+]")[[1]] #)  gsub(" ", "", 
  
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
    f1 <- strsplit(rtermss[grpcheck],":")[[1]]
    f00 <- grep("trait",f1) # position of structure for trait
    if(length(f00)>0){
      ZKgrouping.str <- f1[f00]
    }else{
      ZKgrouping.str <- rep("diag(trait)",length(f1)) 
    }
    f0 <- f1[f00] # structure for eigend
    f2 <- f1[setdiff(1:length(f1),f00)]
    #if(length(f0)>0){f2 <- paste(f0,f2,sep=":")}
    random <- as.formula(paste("~",paste(rtermss[-grpcheck], collapse = " + "))) # new random
    random <- as.formula(paste(as.character(random),collapse=""))
    namere <- expi2(f2) # name of the random effect to look in to grouping argument
    ZKgrouping <- list()
    cous <- 0
    cous2 <- numeric()
    #ZKgrouping.str <- character()
    for(u in 1:length(namere)){
      cous <- cous+1
      cous2[u] <- cous
      grpu <- which(names(grouping)==namere[u])
      Zgrpu <- grouping[[grpu]]
      if(is.null(Zgrpu)){stop("Random effect specified with the grp() function not specified in the grouping argument.\n",call. = FALSE)}
      Gu <- which(names(G)==namere[u])
      if(length(Gu) > 0){
        Kgrp <- G[[Gu]]
        cat("Ignore warning messages about variance-covariance matrices specified in the \nG argument not used for your grouping effects.\n")
      }else{ Kgrp <- diag(ncol(Zgrpu))}
      ZKgrouping[[u]] <- list(Z=Zgrpu, K=Kgrp)
      names(ZKgrouping)[u] <- namere[u]
    }
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
    #i=2
    
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
          j3 <- gsub(".*:","",vara) # part removed
          #j2 <- gsub(" ","", j2)
          k1 <- gsub(":.*","",vara) # to remove in next step
          #regexpr("\\((.*)\\)", j2[1])
          orx <- gsub(k1,"",j2, fixed = TRUE) # order of columns by location
          where <- as.matrix(apply(data.frame(unique(orx)),1,function(x,y){which(y==x)},y=orx)) # each column says indeces for each level of at()
          colnames(where) <- unique(j2)
          for(u in 1:dim(where)[2]){ # u=1
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
    
    if(length(grpcheck)>0){
      #added <- length(Z):(length(Z)+length(ZKgrouping))
      #Z[[added]] <- ZKgrouping
      Z <- c(Z,ZKgrouping)
      cous2 <- max(counter) + cous2
      counterl <- c(counterl,cous2)
      ttt$ran.trt.str <- c(ttt$ran.trt.str,ZKgrouping.str)
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
    
    
    res <- mmer(Y=yvar, X=X, Z=Z, R=R, method=method, init=init,
                iters=iters,tolpar=tolpar,
                tolparinv = tolparinv,draw=draw,silent=silent, 
                constraint = constraint,EIGEND = EIGEND,
                forced=forced,IMP=IMP,complete=complete,restrained=torestrain)
    
    
  }else{###only fixed effects
    res <- glm(yvars~X)
  }
  #########
  #print(str(Z))
  return(res)
}
