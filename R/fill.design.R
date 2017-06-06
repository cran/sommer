fill.design <- function(x,rows="ROW",ranges="RANGE",by){
  
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
  # if user has multiple field
  if(!missing(by)){
    y <- split(x,x[,by])
    xnew <- lapply(y,function(x){
      
      bybo <- na.omit(unique(x[,by]))
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
      return(newf)
    })
    xnew <- do.call(rbind(xnew))
  }else{
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
    xnew <- newf
  }
  return(xnew)
}