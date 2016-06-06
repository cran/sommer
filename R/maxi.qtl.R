maxi.qtl <- function(sma, distan=10, q95=2.5){
  sma <- data.frame(sma)
  param <- vector("list", length(unique(sma$chr)))
  param <- lapply(param, function(x){x <- c(q95,35)})
  ## sma is the dataframe
  ## param is a list with 2 values, lod threshold and #of qtls
  sma <- apply(sma, 2, FUN=as.numeric)
  sma <- as.data.frame(sma)
  if(is.null(param)){
    param <- vector("list", length(unique(sma$chr)))
    param <- lapply(param, function(x){x <- c(1,1)})
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
        to.elim <- distan#abs(prov2$pos[w[k]] - 10)
        zz <- (which(prov2$pos > (provi - to.elim) & prov2$pos < (provi + to.elim) ))[-w[k]]
        if(length(zz) > 0){
          prov2 <- prov2[-zz,]
        }else{prov2 <- prov2}
        selected <- which(prov2$lod == mm)[1]
        prov2 <- prov2[-selected,]
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
  
  return(res2)
  
}