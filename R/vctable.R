vctable <- function(ts,nz,nr){
  traitm <- expand.grid(1:ts,1:ts)
  if(ts > 1){
    traitm <- (traitm[!duplicated(t(apply(traitm, 1, sort))),])[,c(2,1)]; colnames(traitm) <- c("t1","t2")
  }else{
    traitm <- as.matrix(cbind(traitm,traitm)[,1:2]) # when single trait we have issues
    colnames(traitm) <- c("t1","t2")
  }
  
  all <- rep(list(traitm), nz+nr)
  res <- rep(1, nrow(traitm))
  for(u in 2:length(all)){
    all[[u]] <- all[[u]] + max(all[[u-1]])
    res <- c(res,rep(u,nrow(traitm)))
  }
  all2 <- do.call(rbind,all); rownames(all2) <- NULL
  all3 <- cbind(all2,res) # trait combination (t1,t2) for all random effects (res)
  return(all3)
}



vctable.help <- function(random=NULL,rcov=NULL, data){
  
  yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]]
  ttt0<-  apply(data.frame(yuyu),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  
  #ttt0<- strsplit(as.character(random[2]), split = "[+]")[[1]] #) gsub(" ", "", 
  newrandom <- character()
  traitstr <- character()
  for(k in 1:length(ttt0)){ # for each random effect check 
    ttt0.1 <- grep("\\(trait",ttt0[k]) # if user specified the structure for trait
    if(length(ttt0.1)>0){ # if so
      base <- strsplit(ttt0[k],":")[[1]]
      fortrait <- grep("\\(trait",base)
#       res <- setdiff(1:length(base),fortrait) # random effects
#       if(length(res)==2){ # if there was an interaction random effect
#         real.re <- expi(base[res[1]]) # real random effect nae
#       }
      traitstr[k] <- base[fortrait] # the structure for trait
      newrandom[k] <- paste(base[-fortrait],collapse = ":")
    }else{ # if user don't specify the trait structure
      traitstr[k] <- "diag(trait)" # the structure for trait
      newrandom[k] <- ttt0[k]
    }
  }
  random <- as.formula(paste("~",paste(newrandom, collapse = " + "))) # new random
  
  if(is.null(rcov)){
    rcov <- as.formula("~units")
  }
  
  yuyu <- strsplit(as.character(rcov[2]), split = "[+]")[[1]]
  ttt0<-  apply(data.frame(yuyu),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  
  #ttt0<- strsplit(as.character(rcov[2]), split = "[+]")[[1]] #) gsub(" ", "", 
  newrcov <- character()
  traitstr.rcov <- character()
  for(k in 1:length(ttt0)){ # for each random effect check 
    ttt0.1 <- grep("\\(trait",ttt0[k]) # if user specified the structure for trait
    if(length(ttt0.1)>0){ # if so
      base <- strsplit(ttt0[k],":")[[1]]
      fortrait <- grep("\\(trait",base)
      traitstr.rcov[k] <- base[fortrait] # the structure for trait
      newrcov[k] <- paste(base[-fortrait],collapse = ":")
    }else{ # if user don't specify the trait structure
      traitstr.rcov[k] <- "diag(trait)" # the structure for trait
      newrcov[k] <- ttt0[k]
    }
  }
  rcov <- as.formula(paste("~",paste(newrcov, collapse = " + "))) # new random
  res <- list(random=random, rcov=rcov, ran.trt.str=traitstr, rcov.trt.str=traitstr.rcov)
  return(res)
}