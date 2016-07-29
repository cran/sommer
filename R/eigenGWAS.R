eigenGWAS <- function(markers, eivec=1, map=NULL){
  
  ################
  sapply_pb <- function(X, FUN, ...)
  {
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    
    wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir=env)
      setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
    }
    res <- sapply(X, wrapper, ...)
    close(pb)
    res
  }
  
  ################
  apply_pb <- function(X, MARGIN, FUN, ...)
  {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTxtProgressBar(get("pb", envir= env),
                        curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  #################
  if(!is.null(map)){
    required.names <- c("Chrom","Position")
    
    che <- which(names(map)%in%required.names)
    if(length(che) < 2){
      stop("Column names; 'Chrom' and 'Position' need 
           to be present in the data frame provided.",call. = FALSE)
    }
  }
  
  num.check <- which(!is.numeric(markers))
  if(length(num.check) > 0){
    stop("Marker matrix needs to be in numeric format. Please,
         make sure you do not have extra columns or the wrong format.",call. = FALSE)
  }
  NM <- dim(markers)
  N <- NM[1]
  M <- NM[2]
  Cmat=cor(t(markers)) # N X N matrix
  #Cmat<- A.mat(markers)
  Eg=eigen(Cmat)
  COL=c(rep("red", N/2), rep("blue", N/2))
  layout(matrix(c(1,2,2,3,4,4),2,3, byrow=TRUE))
  plot(Eg$vectors[,1], Eg$vectors[,2], col=COL, xlab="Eigen vector 1", ylab="Eigen vector 2", bty='n')
  #EigenGWAS for eigen vector 1
  EV1=scale(Eg$vectors[,1])
  
  #Egwas1=matrix(0, M, 4)
  
  #for(i in 1:M)
  #{
  #  mod=lm(EV1~markers[,i])
  #  Egwas1[i,1]=summary(mod)$coefficients[2,1]
  #  Egwas1[i,2]=summary(mod)$coefficients[2,2]
  #  Egwas1[i,3]=summary(mod)$coefficients[2,4]
  #  Egwas1[i,4]=summary(mod)$coefficients[2,3]^2
  #}
  
  
  Egwas1 <- apply_pb(markers,2,function(x,y){
    
    mod=lm(y~x)
    res <- summary(mod)$coefficients[2,c(1,2,4,3)]
    res[4] <- res[4]^2
    return(res)

  }, y=EV1)
  #str(Egwas1)
  Egwas1 <- t(Egwas1)
  #head(Egwas1)
  #Manhattan plot
  #layout(matrix(1,1,1))
  ########################################################
  # plot the p.values for the EigenGWAS for eigenvector 1
  ########################################################
  #plot((Egwas1[,3]))
  xxx <-matrix(Egwas1[,3])
  rownames(xxx) <- colnames(markers)
  
  
  
  #head(mappo)
  if(!is.null(map)){
    map <- map[which(!duplicated(map$Locus)),]
    rownames(map) <- map$Locus
    mappo <- map[rownames(xxx),]
    mappo$p.val <- NA
    mappo$p.val <- -log10(Egwas1[,3])
    mappo <- mappo[which(!is.na(mappo$Chrom)),]
    
    manhattan(mappo, show.fdr = FALSE)
    title(paste("Signature of selection across loci\nfor eigen vector",eivec))
  }else{
    plot(-log10(Egwas1[,3]), xlab="Loci", ylab=expression(-log[10](italic(p))), bty="n", col="cadetblue", main=paste("Signature of selection across loci\nfor eigen vector",eivec))
  }
  
  
  
  #highlight selected loci
  plot(sort(rchisq(M, 1)), sort(Egwas1[,4]), bty="n", xlab="Expected Chi-sq", ylab="Observed Chi-sq")
  abline(a=0,b=1, col="red")
  
  #fst
  gp=which(EV1>0)
  freqM=colMeans(markers)/2
  freq_G1=colMeans(markers[gp,])/2
  freq_G2=colMeans(markers[-gp,])/2
  w1=length(gp)/N
  w2=1-w1
  Fst=(2*w1*(freq_G1-freqM)^2+w2*(freq_G2-freqM)^2)/(freqM*(1-freqM))
  plot(Fst*N, Egwas1[,4], xlab="Fst", ylab="Chi-sq", bty="n")
  abline(a=0, b=1, col="red")
  
  layout(matrix(1,1,1))
  result <- list(Egwas=-log10(Egwas1[,3]), Fst=Fst)
  return(result)
}





