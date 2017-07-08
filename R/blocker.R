blocker <- function(x,rows="ROWS",ranges="RANGES",by, char=FALSE){
  
  if(missing(by)){
    x$FIELDINST <- "FIELD.1"
    by="FIELDINST"
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
  ### make sure you create provisional rows and cols so the blocking goes well
  
  # if(!missing(by)){
  #   checo3 <- is.unsorted(x[,by])
  #   if(checo){
  #     stop("The dataset provided needs to be sorted by the same argument that you provide in 'by' ",call. = FALSE)
  #   }
  # }
  
  topox <- lapply(topox,function(xx){
    # ranges
    ho1 <- table(xx[,c(rows,ranges)])
    ho2<- 1:length(as.numeric(colnames(ho1)))
    xx$RANGEX <- ho2[as.factor(xx[,ranges])]
    # rows
    ho3 <- 1:length(as.numeric(rownames(ho1)))
    xx$ROWX <- ho3[as.factor(xx[,rows])]
    return(xx)
  })
  
  rows0 <- rows
  ranges0 <- ranges
  
  ranges <- "RANGEX"
  rows <- "ROWX"
  
  lapply(topox,dim)
  names(topox)
  
  rock <- lapply(topox, function(x){ # x <- topox[[5]]
    # potential divisions
    divo <- 2^(seq(1,100))
    divo[1:4]
    
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
    
    #print(divo)
    #print((x[,short]))
    
    divo <- divo[which(divo <= max(x[,short]))]
    
    
    if(length(divo) > 0){
      for(u in 1:length(divo)){
        (useit <- divo[u]) # total number of blocks pursued
        (useit2 <- useit/2)  # how many in row and range direction
        
        if(u ==1){
          useit2 <- useit
        }
        inld <- seq(1,max(x[,large]),floor(max(x[,large]/useit2)))[1:useit2]
        insd <- seq(1,max(x[,short]),floor(max(x[,short]/useit2)))[1:useit2]
        
        inld;insd
        
        inld2 <- (c(inld-1,max(x[,large])))[-1] # final part of the blocks
        insd2 <- (c(insd-1,max(x[,short])))[-1] #
        
        if(u ==1){ # only when 2 blocks are needed, one direction is used
          upo <- numeric(length = dim(x)[1])
          for(o in 1:length(inld)){
            ww <- which(x[,large] >= inld[o] & x[,large] <= inld2[o])
            upo[ww] <- o
          }
          x[,paste("block",useit,sep=".")]<- paste(upo,sep="") 
        }else{ # otherwise we use both directions
          upo <- numeric(length = dim(x)[1])
          for(o in 1:length(inld)){
            ww <- which(x[,large] >= inld[o] & x[,large] <= inld2[o])
            upo[ww] <- o
          }
          
          upo2 <- numeric(length = dim(x)[1])
          for(o in 1:length(insd)){
            ww <- which(x[,short] >= insd[o] & x[,short] <= insd2[o])
            upo2[ww] <- o
          }
          
          ado <- paste(upo,upo2,sep="")
          x[,paste("block",length(unique(ado)),sep=".")]<-  ado
        }
        
      }
      
      
      bolos <-grep("block.",colnames(x))
      bads <- which(apply(as.data.frame(x[,bolos]),2,function(x){length(which(is.na(x)))/length(x)}) == 1)
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
    }
    
    return(x)
  })
  
  ## change to common names among everyone
  rock <- lapply(rock,function(x){
    where <- grep("block.",colnames(x))
    if(length(where)>0){
    colnames(x)[where] <- paste("block.",1:length(where),sep="")
    }
    return(x)
  })
  
  # lapply(rock,function(x){
  #   chocho <- which(colnames(x) %in% "block.2")
  #   if(length(chocho)>0){
  #     desplot(block.2~ROW*RANGE|BLOCK, x, out2 = BLOCK, out1 = SUBBLOCK,
  #             out2.gpar=list(col = "black", lwd = 3, lty = 1), xlab="ROW",ylab="RANGE",
  #             out1.gpar=list(col = "aquamarine1", lwd = 1, lty = 1))
  #   }
  #   
  # })
  
  #library(plyr)
 
  ### 
  ## move back to a good dataframe
  if(!missing(by)){
    
    # common <- lapply(rock,function(x){names(x)}) # names for each data frame
    # common.n <- unique(unlist(common))#Reduce(unique, common) # common names
    # ado3 <- common.n#[grep("block.",common.n)]
    # rock2 <- lapply(rock,function(x,y){ # x <- rock[[2]]
    #   cococo <- setdiff(y,colnames(x)) # we have to add
    #   if(length(cococo)>0){
    #     found <- grep("block.",colnames(x))
    #     if(length(found)>0){
    #       uuu <- colnames(x)[found]
    #       uuu2 <- uuu[length(uuu)]
    #       x[,cococo] <- x[,uuu2]
    #     }else{
    #       x[,cococo] <- NA
    #     }
    #     
    #   }
    #   return(x)
    # },y=ado3) ## add the columns so everyone has the same columns
    #rock <- lapply(rock,function(x,y){x[,y]},y=common.n)
    fin <- do.call(plyr::rbind.fill, rock) # do.call(rbind,rock2) 
  }else{ # if a single field
    fin <- rock[[1]]
  }
  ## remove bad names
 
  ### sort both frames equally
  dim(fin)
  #x <- x[with(x, order(x[,field], x[,rows], x[,ranges])), ]
  # x <- x[ order(x[,field], x[,rows0], x[,ranges0]), ]
  # fin <- fin[ order(fin[,field], fin[,rows], fin[,ranges]), ]
  
  #sum(unlist(lapply(rock,function(x){dim(x)[1]})))
  
  fin <- fin[,setdiff(colnames(fin),c("ROWX","RANGEX"))]
  
  bolos <- colnames(fin)[grep("block.",colnames(fin))]
  
  fin <- merge(x,fin[,c(by,rows0,ranges0,bolos)], by=c(by,rows0,ranges0), all=TRUE)
  
  ## fill NA gaps
  bobo <- grep("block.",colnames(fin))
  if(length(bobo>0)){
    for(tt in 2:length(bobo)){
      missed <- which(is.na(fin[,bobo[tt]]))
      if(length(missed)>0){
        fin[missed,bobo[tt]] <- fin[missed,bobo[tt-1]]
      }
    }
  }
  
  good <- setdiff(colnames(fin),by)
  fin <- fin[,good]
  
  return(fin)
}
