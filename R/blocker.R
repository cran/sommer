blocker <- function(x,rows="ROWS",ranges="RANGES",by, char=FALSE){
  
  checo <- which(colnames(x) %in% c(rows,ranges))
  if(length(checo)!=2){
    stop("Please double check the rows and ranges argument that you provided. We did not find such columns.\n")
  }
  
  roro <- x[,which(colnames(x)==rows)]
  raro <- x[,which(colnames(x)==ranges)]
  
  if(!is.numeric(roro) | !is.numeric(raro)){
    stop("Please make sure that the columns; ",rows," and ",ranges," are both numeric.\n",call. = FALSE)
  }
  
  
  fac <- which(unlist(lapply(x,class))=="factor")
  if(length(fac)>0){
    for(u in 1:length(fac)){
      fac2 <- fac[u]
      x[,fac2] <- as.character(x[,fac2])
    }
  }
  
  if(!missing(by)){
    topox <- split(x,x[,by])
  }else{
    topox <- list(x)
  }
  
  rock <- lapply(topox, function(x){
    # potential divisions
    divo <- 2^(seq(1,100))
    # odd indexes only divide the largest side
    # even indexes diveide both
    
    max(x[,rows])
    max(x[,ranges])
    
    if(max(x[,rows])!=max(x[,ranges])){
      pos <- c(rows,ranges)
      large <- which(c(max(x[,rows]),max(x[,ranges])) == max(c(max(x[,rows]),max(x[,ranges]))))
      short <- setdiff(1:2,large)
      large <- pos[large]
      short <- pos[short]
    }else{
      large <- rows
      short <- ranges
    }
    short
    large
    
    namos <- paste("block",divo,sep=".")
    
    ## just this ones are smaller and therefore to be used for division
    todo <- which(divo < (max(x[,rows])*max(x[,ranges])))
    yeah <- max(todo[-length(todo)])
    yeah <- yeah-1
    
    divo1 <- sort(rep(2^(seq(1,100)),2))
    #divo2 <- sort(c(divo1,divo1-1))
    
    for(l in 1:yeah){
      x[,namos[l]]<- NA
      
      if(l%%2 !=0){ # 2^1=2, 2^3=8
        
        use <- divo1[l]   
        #
        rest <- max(x[,large])%%use
        blocksize <- floor(max(x[,large])/use)
        st <- seq(1,max(x[,large]),blocksize) 
        en <- st+(blocksize-1)
        en[length(en)] <- en[length(en)] + rest # add the rest of rows/ranges that were left
        ## now add block information to all columns
        for(o in 1:length(en)){
          ww <- which(x[,large] >= st[o] & x[,large] <= en[o])
          x[ww,namos[l]] <- o
        }
        
        if(l>1){
          x[,namos[l]] <- paste(x[,namos[l]],x[,namos[l-1]],sep="")
        }
        
      }else if(l%%2 == 0){ # 2^2=4,2^4=16,... 
        
        use <- divo1[l]#/2   
        #
        rest <- max(x[,short])%%use
        blocksize <- floor(max(x[,short])/use)
        st <- (seq(1,max(x[,short]),blocksize))[1:use] 
        en <- (st+(blocksize-1))[1:use]
        en[length(en)] <- en[length(en)] + rest # add the rest of rows/ranges that were left
        ## now add block information to all columns
        for(o in 1:length(en)){
          ww <- which(x[,short] >= st[o] & x[,short] <= en[o])
          x[ww,namos[l]] <- o
        }
        x[,namos[l]] <- paste(x[,namos[l]],x[,namos[l-1]],sep="")
        
      }
      
    }
    
    bolos <-grep("block.",colnames(x))
    bads <- which(apply(x[,bolos],2,function(x){length(which(is.na(x)))/length(x)}) == 1)
    if(length(bads)>0){
      x <- x[-c(bolos[bads])]
    }
    head(x)
    
    bolos2 <-grep("block.",colnames(x))
    bolos2
    for(u in 1:length(bolos2)){
      jj <- bolos2[u]
      x[,jj] <- c(1:10000)[as.factor(as.character(x[,jj]))]
    }
    head(x)
    cat("Number of levels achieved.\n")
    print(apply(x[,bolos],2,function(x){length(unique(x))}))
    ## return a plot that gives them an idea
    #layout(matrix(1:6,2,3))
#     for(u in 1:length(bolos2)){
#       jj <- bolos2[u]
#       pp <- levelplot(as.formula(paste(names(x)[jj],"~",rows,"*",ranges)),main=names(x)[jj], data=x)
#       print(pp)
#     }
    ## if user wants them as character
    if(char){
      for(u in 1:length(bolos2)){
        jj <- bolos2[u]
        x[,jj] <- paste("B",x[,jj],sep=".")
      }
    }
    return(x)
  })
  ### 
  ## move back to a good dataframe
  if(!missing(by)){
    fin <- do.call(rbind,rock)
  }else{
    fin <- rock[[1]]
  }
  return(fin)
}