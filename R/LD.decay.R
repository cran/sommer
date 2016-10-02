LD.decay <- function(markers,map,silent=FALSE,unlinked=FALSE,gamma=.95){
  
  #if(is.null(rownames(map))){
  good <- which(!duplicated(map$Locus))
  map <- map[good,]
  rownames(map) <- map$Locus 
  #}
  ## clean markers and map from sd=0 and markers with no position
  markers <- markers[,which(apply(markers,2,sd)>0)]
  map <- map[which(!is.na(map$LG)),] # which(apply(map,1,function(x){length(which(is.na(x)))}) ==0)
  #############
  fullmap <- map
  lgs <- unique(map$LG)
  ############################
  dat.list <- list(NA) # will contain the LD (r) and genetic distance (d) for each LG
  ############################
  if (!silent) {
    count <- 0
    tot <- length(lgs)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  if(!unlinked){
    
    LDM <- list()
    for(k in lgs){ # for each linkaghe group
      if (!silent) {
        count <- count + 1
      }
      ords2 <- fullmap[which(fullmap[,"LG"] == k),]
      #ords2 <- ords
      head(ords2);tail(ords2)
      #################
      #intersect markers and map
      #################  
      inn <- intersect(ords2$Locus,colnames(markers))
      
      if(length(inn)>0){
        ords2 <- ords2[inn,]
        
        cor2.mat <- cor(as.matrix(markers[,inn]),use="pairwise.complete.obs")^2
        LDM[[k]] <- cor2.mat
        
        N <- dim(markers[,inn])[1]
        cor2.matp <- 1-pchisq(cor2.mat*N,df=1)
        
        #cor(fofo[,which(apply(fofo,2,sd) > 0)])^2
        # trnsforming to double haploid data
        # this is a mtrix of genetic distances
        mat.dist <- matrix(0,dim(ords2)[1], dim(ords2)[1])
        for(i in 1:dim(ords2)[1]){
          for(j in 1:dim(ords2)[1]){
            mat.dist[i,j] <- ords2$Position[i] - ords2$Position[j]
          }
        }
        # get absolute values
        mat.dist <- abs(mat.dist) 
        # make zero the valuies of upper and lower triangular
        cor2.mat[lower.tri(cor2.mat)] <- 0
        cor2.matp[lower.tri(cor2.mat)] <- 0
        mat.dist[lower.tri(mat.dist)] <- 0
        
        y.dist <- as.vector(mat.dist) #los deshace columna por column
        y.cor <- as.vector(cor2.mat) 
        p.cor <- as.vector(cor2.matp)
        
        dat <- data.frame(d=y.dist,r2=y.cor,p=p.cor)
        head(dat)
        dat2 <- dat[-which(dat$d == 0 & dat$r == 0),]
        dat.list[[k]] <- dat2
        if (!silent) {
          setTxtProgressBar(pb, (count/tot))
        }
      }else{
        dat <- data.frame(d=0,r2=0,p=0)
        head(dat)
        dat2 <- dat[-which(dat$d == 0 & dat$r == 0),]
        dat.list[[k]] <- dat2
      }
      
    }
    # make a big data frame
    if (!silent) {
      setTxtProgressBar(pb, (count/tot))
    }
    
    big <- do.call("rbind",dat.list)
    
    # big contains all the LD and and genetic distances for all groups
    resp <- list(by.LG=dat.list, all.LG=big, LDM=LDM)
    
  }else{ # if WANTS TO KNOW THE THRESHOLD OF UNLINKED
    doo <- intersect(colnames(markers), fullmap$Locus)
    cor2.mat <- cor(as.matrix(markers[,doo]),use="pairwise.complete.obs")^2
    
    dat.list <- numeric()#list()#store the r2 values of unlinked markers with chromosome k
    for(k in lgs){ # for each linkaghe group
      if (!silent) {
        count <- count + 1
      }
      # markers in kth chromosome
      ords2 <- fullmap[which(fullmap[,"LG"] == k),]
      # markers not in kth chromosome
      ords3 <- fullmap[which(fullmap[,"LG"] != k),]
      # markers in kth and present
      sik <- intersect(ords2$Locus,colnames(cor2.mat))
      # markers not in kth and present
      nok <- intersect(ords3$Locus,colnames(cor2.mat))
      #ords2 <- ords
      #head(ords2);tail(ords2)
      #into <- setdiff(nok,ords2$Locus)
      
      step1 <- cor2.mat[sik,nok]
      #boxplot(unlist(step1))
      dat.list[k] <- quantile(sqrt(step1[upper.tri(step1)]),gamma)^2
      
      if (!silent) {
        setTxtProgressBar(pb, (count/tot))
      }
    }
    # make a big data frame
    resp <- list(by.LG=dat.list, all.LG=mean(dat.list))
  }
  return(resp)
}