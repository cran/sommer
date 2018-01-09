blocker <- function(dat,rows="ROW",ranges="RANGE",by=NULL){
  
  if(is.null(by)){by="FIELDINST"; dat[,"FIELDINST"] <- "F1"}
  
  checo <- which(colnames(dat) %in% c(rows,ranges))
  if(length(checo) != 2){
    stop("Please double check the rows and ranges argument that you provided. We did not find such columns.\n")
  }
  
  roro <- dat[,which(colnames(dat)==rows)]
  raro <- dat[,which(colnames(dat)==ranges)]
  
  if(!is.numeric(roro) | !is.numeric(raro)){
    stop("Please make sure that the columns; ",rows," and ",ranges," are both numeric.\n",call. = FALSE)
  }
  
  
  fac <- which(unlist(lapply(dat,class))=="factor")
  if(length(fac)>0){
    for(u in 1:length(fac)){
      fac2 <- fac[u]
      dat[,fac2] <- as.character(dat[,fac2])
    }
  }
  
  if(!missing(by)){
    topox <- split(dat,dat[,by])
  }else{
    topox <- list(dat)
  }
  ### make sure you create provisional rows and cols so the blocking goes well
  
  # if(!missing(by)){
  #   checo3 <- is.unsorted(dat[,by])
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
  
  rock <- lapply(topox, function(dat){ # dat <- topox[[5]]
    # potential divisions
    divo <- 2^(seq(1,100))
    divo[1:4]
    
    # odd indexes only divide the largest side
    # even indexes diveide both
    
    max(dat[,rows])
    max(dat[,ranges])
    
    
    if(max(dat[,rows])!=max(dat[,ranges])){
      pos <- c(rows,ranges)
      large <- which(c(max(dat[,rows]),max(dat[,ranges])) == max(c(max(dat[,rows]),max(dat[,ranges]))))
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
    #print((dat[,short]))
    
    divo <- divo[which(divo <= max(dat[,short]))]
    
    
    if(length(divo) > 0){
      for(u in 1:length(divo)){
        (useit <- divo[u]) # total number of blocks pursued
        (useit2 <- useit/2)  # how many in row and range direction
        
        if(u ==1){
          useit2 <- useit
        }
        inld <- seq(1,max(dat[,large]),floor(max(dat[,large]/useit2)))[1:useit2]
        insd <- seq(1,max(dat[,short]),floor(max(dat[,short]/useit2)))[1:useit2]
        
        inld;insd
        
        inld2 <- (c(inld-1,max(dat[,large])))[-1] # final part of the blocks
        insd2 <- (c(insd-1,max(dat[,short])))[-1] #
        
        if(u ==1){ # only when 2 blocks are needed, one direction is used
          upo <- numeric(length = dim(dat)[1])
          for(o in 1:length(inld)){
            ww <- which(dat[,large] >= inld[o] & dat[,large] <= inld2[o])
            upo[ww] <- o
          }
          dat[,paste("block",useit,sep=".")]<- paste(upo,sep="") 
        }else{ # otherwise we use both directions
          upo <- numeric(length = dim(dat)[1])
          for(o in 1:length(inld)){
            ww <- which(dat[,large] >= inld[o] & dat[,large] <= inld2[o])
            upo[ww] <- o
          }
          
          upo2 <- numeric(length = dim(dat)[1])
          for(o in 1:length(insd)){
            ww <- which(dat[,short] >= insd[o] & dat[,short] <= insd2[o])
            upo2[ww] <- o
          }
          
          ado <- paste(upo,upo2,sep="")
          dat[,paste("block",length(unique(ado)),sep=".")]<-  ado
        }
        
      }
      
      
      bolos <-grep("block.",colnames(dat))
      bads <- which(apply(as.data.frame(dat[,bolos]),2,function(dat){length(which(is.na(dat)))/length(dat)}) == 1)
      if(length(bads)>0){
        dat <- dat[-c(bolos[bads])]
      }
      head(dat)
      
      bolos2 <-grep("block.",colnames(dat))
      bolos2
      for(u in 1:length(bolos2)){
        jj <- bolos2[u]
        dat[,jj] <- c(1:10000)[as.factor(as.character(dat[,jj]))]
      }
      head(dat)

      
      ## return a plot that gives them an idea
      #layout(matrix(1:6,2,3))
      #     for(u in 1:length(bolos2)){
      #       jj <- bolos2[u]
      #       pp <- levelplot(as.formula(paste(names(dat)[jj],"~",rows,"*",ranges)),main=names(dat)[jj], data=dat)
      #       print(pp)
      #     }
      ## if user wants them as character
      # if(char){
      #   for(u in 1:length(bolos2)){
      #     jj <- bolos2[u]
      #     dat[,jj] <- paste("B",dat[,jj],sep=".")
      #   }
      # }
    }
    
    return(dat)
  })
  
  ## change to common names among everyone
  rock <- lapply(rock,function(xd){
    where <- grep("block.",colnames(xd))
    if(length(where)>0){
    colnames(xd)[where] <- paste("block.",1:length(where),sep="")
    }
    return(xd)
  })
  
  # lapply(rock,function(dat){
  #   chocho <- which(colnames(dat) %in% "block.2")
  #   if(length(chocho)>0){
  #     desplot(block.2~ROW*RANGE|BLOCK, dat, out2 = BLOCK, out1 = SUBBLOCK,
  #             out2.gpar=list(col = "black", lwd = 3, lty = 1), xlab="ROW",ylab="RANGE",
  #             out1.gpar=list(col = "aquamarine1", lwd = 1, lty = 1))
  #   }
  #   
  # })
  
  #library(plyr)
 
  ### 
  ## move back to a good dataframe
  if(!missing(by)){
    
    # common <- lapply(rock,function(dat){names(dat)}) # names for each data frame
    # common.n <- unique(unlist(common))#Reduce(unique, common) # common names
    # ado3 <- common.n#[grep("block.",common.n)]
    # rock2 <- lapply(rock,function(dat,y){ # dat <- rock[[2]]
    #   cococo <- setdiff(y,colnames(dat)) # we have to add
    #   if(length(cococo)>0){
    #     found <- grep("block.",colnames(dat))
    #     if(length(found)>0){
    #       uuu <- colnames(dat)[found]
    #       uuu2 <- uuu[length(uuu)]
    #       dat[,cococo] <- dat[,uuu2]
    #     }else{
    #       dat[,cococo] <- NA
    #     }
    #     
    #   }
    #   return(dat)
    # },y=ado3) ## add the columns so everyone has the same columns
    #rock <- lapply(rock,function(dat,y){dat[,y]},y=common.n)
    fin <- do.call(plyr::rbind.fill, rock) # do.call(rbind,rock2) 
  }else{ # if a single field
    fin <- rock[[1]]
  }
  ## remove bad names
 
  ### sort both frames equally
  dim(fin)
  #dat <- dat[with(dat, order(dat[,field], dat[,rows], dat[,ranges])), ]
  # dat <- dat[ order(dat[,field], dat[,rows0], dat[,ranges0]), ]
  # fin <- fin[ order(fin[,field], fin[,rows], fin[,ranges]), ]
  
  #sum(unlist(lapply(rock,function(dat){dim(dat)[1]})))
  
  fin <- fin[,setdiff(colnames(fin),c("ROWX","RANGEX"))]
  
  bolos <- colnames(fin)[grep("block.",colnames(fin))]
  
  fin <- merge(dat,fin[,c(by,rows0,ranges0,bolos)], by=c(by,rows0,ranges0), all=TRUE)
  
  ## fill NA gaps
  bobo <- grep("block.",colnames(fin))
  if(length(bobo)>1){
    for(tt in 2:length(bobo)){
      missed <- which(is.na(fin[,bobo[tt]]))
      if(length(missed)>0){
        fin[missed,bobo[tt]] <- fin[missed,bobo[tt-1]]
      }
    }
  }
  
  return(fin)
}

############################


blockerL <- function(dat, nr= 5, rows="ROW",ranges="RANGE", by=NULL,
                     shiftF=0, shiftB=0){
  
  if(is.null(by)){by="FIELDINST"; dat[,"FIELDINST"] <- "F1"}
  
  
  direction <- rows
  ndirection <- ranges
  
  dtlist <- split(dat, dat[,by])
  
  dtlist <- lapply(dtlist, function(dt){
    
    stf <- seq( (1+shiftF) - (nr * 8), max(dt[, rows]), nr) # forward blocking direction
    (heigth <- nr)#abs(round(((max(dt[,rows])*2)*sin(d)*sin(d))/sin((90-d)*2))))
    counter <- 0
    
    ################# forward blocking
    #################
    dt$LBLOCKF <- NA
    
    for(u in 1:length(stf)){ # for each triangule formed
      counter <- counter+1
      (rowsinblock <- stf[u]:max(dt[,rows])) # rows in block
      for(k in 1:length(rowsinblock)){ # for each row in the new block find height and add plot information
        (heigth0 <- max(dt[,ndirection]))#abs(round(((k*2)*sin(d)*sin(d))/sin((90-d)*2))))
        found <- which(dt[,rows] >= rowsinblock[k] & dt[,ranges] <= heigth0)
        if(length(found)>0){dt[found,"LBLOCKF"] <- u}
      }
    }
    #levelplot(TBLOCK~Row*Col, data=dt)
    
    ################# forward blocking
    #################
    dt$LBLOCKB <- NA
    
    stf <- seq( (1+shiftB) - (nr * 8), max(dt[, ranges]), nr) # forward blocking direction
    (heigth <- nr)#abs(round(((max(dt[,rows])*2)*sin(d)*sin(d))/sin((90-d)*2))))
    counter <- 0
    
    for(u in 1:length(stf)){ # for each triangule formed
      counter <- counter+1
      (rowsinblock <- stf[u]:max(dt[,ranges])) # rows in block
      for(k in 1:length(rowsinblock)){ # for each row in the new block find height and add plot information
        (heigth0 <- max(dt[,direction]))#abs(round(((k*2)*sin(d)*sin(d))/sin((90-d)*2))))
        found <- which(dt[,ranges] >= rowsinblock[k] & dt[,rows] <= heigth0)
        if(length(found)>0){dt[found,"LBLOCKB"] <- u}
      }
    }
    
    ################
    # ################ backwards blocking
    # dt$LBLOCKB <- NA
    # (stb <- sort(seq(1,max(dt[,ranges])+(nr*8),nr), decreasing = TRUE)) # backwards
    # for(u in 1:length(stb)){ # for each triangule formed
    #   counter <- counter+1
    #   (rowsinblock <- min(dt[,ranges]):stb[u]) # rows in block
    #   (rowsinblock2 <- stb[u]:min(dt[,ranges])) # to make sure that the smalles height is used at the last rows
    #   for(k in 1:length(rowsinblock)){ # for each row in the new block find height and add plot information
    #     (heigth0 <- max(dat[,direction]))# abs(round(((k*2)*sin(d)*sin(d))/sin((90-d)*2))))
    #     found <- which(dt[,ranges] <= rowsinblock2[k] & dt[,rows] <= heigth0)
    #     if(length(found)>0){dt[found,"LBLOCKB"] <- u}
    #   }
    # }
    
    ## check that there's no blocks smaller than 4 plots
    # nblockf <- unique(dt[,"LBLOCKF"])
    # for(k in nblockf){
    #   kkk <- which(dt[,"LBLOCKF"] == k)
    #   if(length(kkk) < 4){
    #     dt[kkk,"LBLOCKF"] <- dt[kkk,"LBLOCKF"] + 1
    #   }
    #   kkk <- which(dt[,"LBLOCKF"] == k)
    #   if(length(kkk) < 4){
    #     dt[kkk,"LBLOCKF"] <- dt[kkk,"LBLOCKF"] - 2
    #   }
    # }
    # ## check that there's no blocks smaller than 4 plots for backwards
    # nblockb <- unique(dt[,"LBLOCKB"])
    # for(k in nblockb){
    #   kkk <- which(dt[,"LBLOCKB"] == k)
    #   if(length(kkk) < 4){
    #     dt[kkk,"LBLOCKB"] <- dt[kkk,"LBLOCKB"] + 1
    #   }
    #   kkk <- which(dt[,"LBLOCKB"] == k)
    #   if(length(kkk) < 4){
    #     dt[kkk,"LBLOCKB"] <- dt[kkk,"LBLOCKB"] - 2
    #   }
    # }
    return(dt)
  })
  dtt <- do.call(rbind,dtlist)
  #levelplot(TBLOCK~Row*Col, data=dd2)
  return(dtt)
}


#######################







blockerT <- function(dat, d= 60, nr= 5, rows="ROW",ranges="RANGE", 
                     by=NULL, shiftF=0, shiftB=0){
  
  if(is.null(by)){by="FIELDINST"; dat[,"FIELDINST"] <- "F1"}
  
  dtlist <- split(dat, dat[,by])
  
  dtlist <- lapply(dtlist, function(dt){
    
    stf <- seq((shiftF+1) - (nr * 40), max(dt[, rows]), nr) # forward blocking direction
    (heigth <- abs(round(((max(dt[,rows])*2)*sin(d)*sin(d))/sin((90-d)*2))))
    counter <- 0
    
    ################# forward blocking
    #################
    dt$TBLOCKF <- NA
    
    for(u in 1:length(stf)){ # for each triangule formed
      counter <- counter+1
      (rowsinblock <- stf[u]:max(dt[,rows])) # rows in block
      for(k in 1:length(rowsinblock)){ # for each row in the new block find height and add plot information
        (heigth0 <- abs(round(((k*2)*sin(d)*sin(d))/sin((90-d)*2))))
        found <- which(dt[,rows] >= rowsinblock[k] & dt[,ranges] <= heigth0)
        if(length(found)>0){dt[found,"TBLOCKF"] <- u}
      }
    }
    #levelplot(TBLOCK~Row*Col, data=dt)
    
    ################
    ################ backwards blocking
    dt$TBLOCKB <- NA
    (stb <- sort(seq((shiftB+1), max(dt[, rows]) + (nr * 40), nr), ## ###!!!!!!
                 decreasing = TRUE))
    for(u in 1:length(stb)){ # for each triangule formed
      counter <- counter+1
      (rowsinblock <- min(dt[,rows]):stb[u]) # rows in block
      (rowsinblock2 <- stb[u]:min(dt[,rows])) # to make sure that the smalles height is used at the last rows
      for(k in 1:length(rowsinblock)){ # for each row in the new block find height and add plot information
        (heigth0 <- abs(round(((k*2)*sin(d)*sin(d))/sin((90-d)*2))))
        found <- which(dt[,rows] < rowsinblock2[k] & dt[,ranges] < heigth0)
        if(length(found)>0){dt[found,"TBLOCKB"] <- u}
      }
    }
    
    ## check that there's no blocks smaller than 4 plots
    nblockf <- unique(dt[,"TBLOCKF"])
    for(k in nblockf){
      kkk <- which(dt[,"TBLOCKF"] == k)
      if(length(kkk) < 4){
        dt[kkk,"TBLOCKF"] <- dt[kkk,"TBLOCKF"] + 1
      }
      kkk <- which(dt[,"TBLOCKF"] == k)
      if(length(kkk) < 4){
        dt[kkk,"TBLOCKF"] <- dt[kkk,"TBLOCKF"] - 2
      }
    }
    ## check that there's no blocks smaller than 4 plots for backwards
    nblockb <- unique(dt[,"TBLOCKB"])
    for(k in nblockb){
      kkk <- which(dt[,"TBLOCKB"] == k)
      if(length(kkk) < 4){
        dt[kkk,"TBLOCKB"] <- dt[kkk,"TBLOCKB"] + 1
      }
      kkk <- which(dt[,"TBLOCKB"] == k)
      if(length(kkk) < 4){
        dt[kkk,"TBLOCKB"] <- dt[kkk,"TBLOCKB"] - 2
      }
    }
    return(dt)
  })
  dtt <- do.call(rbind,dtlist)
  #levelplot(TBLOCK~Row*Col, data=dd2)
  return(dtt)
}