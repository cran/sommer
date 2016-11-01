maxi.qtl <-   function(sma, distan=10, no.qtl=5, q95=2.5, LOD.int=FALSE, LODdrop=2, allCI=TRUE){
  param=NULL
  sma <- data.frame(sma)
  rnn <- rownames(sma)
  param <- vector("list", length(unique(sma$chr)))
  param <- lapply(param, function(x){x <- c(q95,no.qtl)})
  ## sma is the dataframe
  ## param is a list with 2 values, lod threshold and #of qtls
  sma <- apply(sma, 2, FUN=as.numeric)
  sma <- as.data.frame(sma)
  rownames(sma) <- rnn
  if(is.null(param)){
    param <- vector("list", length(unique(sma$chr)))
    param <- lapply(param, function(x){x <- c(1,1)})
  }
  
  big.peaks.col <- function(x, tre){
    r <- rle(x)
    v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
    pos <- v[which(x[v] > tre)] #position of the real peaks
    hei <- x[pos]
    out <- list(pos=pos,hei=hei)
    return(out)
  }
  
  result <- list(NA)
  for(j in 1:max(sma$chr,na.rm=TRUE)){
    prov <- sma[which(sma$chr == j),]
    para <- param[[j]]
    roxy <- big.peaks.col(prov$lod, tre=para[1]) # first element is the lod threshold, 2nd is the #of qtls
    if(length(roxy$pos) > 0){
      prov2.1 <- prov[roxy$pos,] # to keep
      prov2 <- prov[roxy$pos,] # to discard
      
      w <- numeric()
      res <- data.frame(matrix(NA,nrow=para[2], ncol=3)) # 3 columns because always rqtl gives that
      names(res) <- names(prov2)
      
      GOOD <- TRUE
      k=1
      while(GOOD){
        #for(k in 1:para[2]){# for the peaks required
        mm <- max(prov2$lod, na.rm=TRUE)
        w[k] <- which(prov2$lod == mm)[1]
        provi <- prov2$pos[w[k]]
        
        res[k,c(1:3)] <- prov2[w[k],c(1:3)]
        rownames(res)[k] <- rownames(prov2)[w[k]] #$$$$$$$$$$$$$$$$$$$$
        to.elim <- distan#abs(prov2$pos[w[k]] - 10)
        zz <- (which(prov2$pos > (provi - to.elim) & prov2$pos < (provi + to.elim) ))[-w[k]]
        if(length(zz) > 0){
          prov2 <- prov2[-zz,]
        }else{prov2 <- prov2}
        selected <- which(prov2$lod == mm)[1]
        prov2 <- prov2[-selected,]
        #rownames(prov2.1)[-selected]
        k=k+1
        if(length(which(is.na(prov2))) > 0 | dim(prov2)[1] == 0 ){GOOD<-FALSE}
        #}
      }
      
      #maxqtls <- sort(prov2$lod, decreasing = TRUE)
      result[[j]] <- res#prov2[which(prov2$lod %in% maxqtls[1:para[2]]),]
    }
    
  }# end for each chromosome
  #############
  result <- lapply(result, function(x){unique(x)})
  good <- which(unlist(lapply(result, is.null)) == FALSE)
  result <- result[good]
  
  good <- which(unlist(lapply(result, function(m){is.null(dim(m))})) == FALSE)
  result <- result[good]
  
  res2 <- do.call("rbind", result)
  
  res2<- res2[which(!is.na(res2$chr)),]
  
  res3 <- res2
  if(LOD.int){
    
    #print(res3)
    
    lod2 <-list() 
    for(i in 1:dim(res3)[1]){
      #apply(res3,1,function(x,sma){
      
      x <- res3[i,]
      #print(x)
      mpp <- sma[which(sma$chr == as.numeric(x[1])),]
      
      babo <- which(rownames(mpp) == rownames(x))#which(mpp$chr == as.numeric(x[1]) & mpp$pos == as.numeric(x[2]))
      toch <- mpp[babo,]
      baba <- babo-1
      babu <- babo+1
      rt=0
      #print(toch)
      while((rt < LODdrop) & (baba > 1)){ # 2-lod interval
        rt <- abs(as.numeric(mpp[baba,] - toch)[3])
        baba <- baba - 1
      }
      ## baba bow tell us what is the lower bound
      mpp[baba,]
      
      rt=0
      while((rt < LODdrop) & (babu < dim(mpp)[1])){
        rt <- abs(as.numeric(mpp[babu,] - toch)[3])
        babu <- babu + 1
      }
      mpp[babu,]
      if(allCI){
        resx <- rbind( mpp[baba:(babo-1),], toch,  mpp[(babo+1):babu,])
      }else{
        resx <- rbind( mpp[baba,], toch,  mpp[babu,])
      }
      
      
      lod2[[i]] <- resx
    }
    
    res2 <- do.call("rbind",lod2)
  }
  return(res2)
  
}