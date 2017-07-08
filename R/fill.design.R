
fill.design <- function(x,rows="ROW",ranges="RANGE", by, extra){
  
  checo <- which(colnames(x) %in% c(rows,ranges))
  if(length(checo)!=2){
    stop("Please double check the rows and ranges argument that you provided. We did not find such columns.\n")
  }
  
  roro <- x[,which(colnames(x)==rows)]
  raro <- x[,which(colnames(x)==ranges)]
  
  if(!is.numeric(roro) | !is.numeric(raro)){
    stop("Please make sure that the columns; ",rows," and ",ranges," are both numeric.\n",call. = FALSE)
  }
  
  if(!missing(extra)){
    extras=TRUE
  }else{extras=FALSE}
  
  fac <- which(unlist(lapply(x,class))=="factor")
  if(length(fac)>0){
    for(u in 1:length(fac)){
      fac2 <- fac[u]
      x[,fac2] <- as.character(x[,fac2])
    }
  }
  
  
  # if user has multiple field
  if(!missing(by)){ # if multiple fields
    y <- split(x,x[,by])
    xnew <- lapply(y,function(x, extra, extras){ # the data, the additional element, and the T/F
      
      x <- x[order(x[,rows], x[,ranges]), ] # order by rows and columns before starint the whole process
      
      roro <- x[,which(colnames(x)==rows)]
      raro <- x[,which(colnames(x)==ranges)]
      
      bybo <- na.omit(unique(x[,by])) # level of by
      ## rows needed
      ro1 <- seq(min(roro,na.rm=TRUE),max(roro,na.rm=TRUE))
      ## ranges needed
      ra1 <- seq(min(raro,na.rm=TRUE),max(raro,na.rm=TRUE))
      ## order the needed one
      needed <- expand.grid(ro1,ra1); colnames(needed) <- c(rows,ranges)
      needed <- needed[ order(needed[,rows], needed[,ranges]), ]
      head(needed)
      ## order the provided
      x <- x[ order(x[,rows], x[,ranges]), ]
      head(x)
      
      ## create empty frame
      dis <- dim(x)
      dis2 <- dim(needed)
      newf <- data.frame(matrix(NA,dis2[1],dis[2]))
      colnames(newf) <- colnames(x)
      head(newf)
      ## create rownames in both so is easier to merge
      
      #1) store rownames of original dataset
      sto <- rownames(x)
      #2) assign new rownames
      rownames(x) <- paste(x[,rows],x[,ranges],sep=".")
      
      #3) do the same for the new frame
      rownames(needed) <- rownames(newf) <- paste(needed[,rows],needed[,ranges],sep=".")
      
      #4) complete the merge
      newf[rownames(x),] <- x
      newf[,c(rows,ranges)] <- needed
      head(newf)
      rownames(newf) <- NULL
      newf[,by] <- bybo
      
      
      with(newf, table(FIELDINST, BLOCK, useNA = "always"))
      dim(newf)
      with(newf, table(FIELDINST, ROW, useNA = "always"))
      with(newf, table(FIELDINST, RANGE, useNA = "always"))
      
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      if(extras){
        for(u in 1:length(extra)){
          pextra <- extra[u]
          pox <- table(newf[,c(rows,ranges,pextra)]) # all zero should be filled in
          leve <- dim(pox)[3] # levels for the extra factor
          levelnames <- na.omit(unique(newf[,pextra]))
          ################
          for(o in 1:leve){ # for each level of the xtra column figure out which rows and ranges are there
            ## let's see until which row extends
            init1 <- which(apply(pox[,,o],1,sum) > 0) #end in row direction
            stend1 <- c(init1[1],init1[length(init1)]) # start and end in row direction
            
            init2 <- which(apply(pox[,,o],2,sum) > 0) #end in row direction
            stend2 <- c(init2[1],init2[length(init2)]) # start and end in range direction
            
            kk <- which(newf[,rows] >= stend1[1] & newf[,rows] <= stend1[2] & newf[,ranges] >= stend2[1] & newf[,ranges] <= stend2[2] )
            newf[kk,pextra] <- levelnames[o]
          }
          ################
        }
      }
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      return(newf)
    }, extra=extra, extras=extras)
    
    xnew <- do.call(rbind,xnew)
  }else{
    
    roro <- x[,which(colnames(x)==rows)]
    raro <- x[,which(colnames(x)==ranges)]
    
    cat("Argument 'by' not provided. Single field assumed.\n")
    ro1 <- seq(min(roro,na.rm=TRUE),max(roro,na.rm=TRUE))
    ## ranges needed
    ra1 <- seq(min(raro,na.rm=TRUE),max(raro,na.rm=TRUE))
    ## order the needed one
    needed <- expand.grid(ro1,ra1); colnames(needed) <- c(rows,ranges)
    needed <- needed[ order(needed[,rows], needed[,ranges]), ]
    head(needed)
    ## order the provided
    x <- x[ order(x[,rows], x[,ranges]), ]
    head(x)
    
    ## create empty frame
    dis <- dim(x)
    dis2 <- dim(needed)
    newf <- data.frame(matrix(NA,dis2[1],dis[2]))
    colnames(newf) <- colnames(x)
    head(newf)
    ## create rownames in both so is easier to merge
    
    #1) store rownames of original dataset
    sto <- rownames(x)
    #2) assign new rownames
    rownames(x) <- paste(x[,rows],x[,ranges],sep=".")
    
    #3) do the same for the new frame
    rownames(needed) <- rownames(newf) <- paste(needed[,rows],needed[,ranges],sep=".")
    
    #4) complete the merge
    newf[rownames(x),] <- x
    newf[,c(rows,ranges)] <- needed
    head(newf)
    rownames(newf) <- NULL
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    if(extras){
      for(u in 1:length(extra)){
        pextra <- extra[u]
        pox <- table(newf[,c(rows,ranges,pextra)]) # all zero should be filled in
        leve <- dim(pox)[3] # levels for the extra factor
        levelnames <- na.omit(unique(newf[,pextra]))
        ################
        for(o in 1:leve){ # for each level of the xtra column figure out which rows and ranges are there
          ## let's see until which row extends
          init1 <- which(apply(pox[,,o],1,sum) > 0) #end in row direction
          stend1 <- c(init1[1],init1[length(init1)]) # start and end in row direction
          
          init2 <- which(apply(pox[,,o],2,sum) > 0) #end in row direction
          stend2 <- c(init2[1],init2[length(init2)]) # start and end in range direction
          
          kk <- which(newf[,rows] >= stend1[1] & newf[,rows] <= stend1[2] & newf[,ranges] >= stend2[1] & newf[,ranges] <= stend2[2] )
          newf[kk,pextra] <- levelnames[o]
        }
        ################
      }
    }
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    
    xnew <- newf
  }
  
  return(xnew)
}


