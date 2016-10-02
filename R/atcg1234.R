atcg1234 <- function(data, ploidy=2, format="ATCG", maf=0, multi=TRUE, silent=FALSE, by.allele=FALSE){
  ##### START GBS.TO.BISNP DATA ######
  gbs.to.bisnp <- function(x) {
    y <- rep(NA,length(x))
    y[which(x=="A")] <- "AA"
    y[which(x=="T")] <- "TT"
    y[which(x=="C")] <- "CC"
    y[which(x=="G")] <- "GG"
    y[which(x=="R")] <- "AG"
    y[which(x=="Y")] <- "CT"
    y[which(x=="S")] <- "CG"
    y[which(x=="W")] <- "AT"
    y[which(x=="K")] <- "GT"
    y[which(x=="M")] <- "AC"
    y[which(x=="+")] <- "++"
    y[which(x=="0")] <- "NN"
    y[which(x=="-")] <- "--"
    y[which(x=="N")] <- NA
    return(y)
  }
  ##### END GBS.TO.BISNP DATA ######
  imputeSNP <- function(data){
    #######
    data2 <- apply(data,2,function(x){
      areNA <- which(is.na(x))
      if(length(areNA)>0){
        pos.all <- table(data[,1])
        totake <- names(pos.all)[which(pos.all == max(pos.all))]
        x[areNA] <- totake
      }
      return(x)
    })
    #######
    return(data2)
  }
  #### apply with progress bar ######
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
  ###### zero.one function
  zero.one <- function(da){
    # this function takes a matrix of markers in biallelic format and returns a matrix of
    # presense/absense of alleles
    mar.nam <- colnames(da)#unique(gsub("\\.\\d","", names(da))) # find a dot and a number after the dot
    mat.list <- list(NA) # list of matrices for each marker
    wi=0 # counter
    if(!silent){
      count <- 0
      tot <- length(mar.nam)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    for(i in 1:length(mar.nam)){ # for each marker
      wi=wi+1
      if(!silent){
        count <- count + 1
      }
      
      v <- which(colnames(da)==mar.nam[i])#grep(mar.nam[i], colnames(da))
      
      if(length(v)==0){
        qqqqq <- grep(mar.nam[i-1],names(da))
        qqqqq2 <- names(da)[qqqqq[length(qqqqq)] + 1]
        
        stop(paste("Marker",qqqqq2,"has a problem"), call.=FALSE)
      }else if(length(v) == 1){ # for markers with a single column
        prov <- matrix(da[,v])
      }else{prov <- da[,v]}
      ##################################
      alls <- unique(unlist(strsplit(prov,"")))
      alls <- alls[which(!is.na(alls))]
      ninds <- dim(prov)[1]
      fff <- apply(data.frame(alls),1,function(h){
        temp <- numeric(length = ninds)
        temp[grep(h,prov)]<-1
        #make sure is full rank
        
        return(temp)
      })#1 # assigning 1's
      #if(FULL){ # if user want to make sure only get the columns that will ensure full rank
      #  fff <- t(unique(t(fff)))
      #}
      colnames(fff) <- paste(mar.nam[i],alls, sep="/")
      
      mat.list[[i]] <- fff
      if(!silent){
        setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
      }
      
    }
    
    fin.mat <- do.call(cbind,mat.list)
    rownames(fin.mat) <- rownames(da)
    
    # make full rank
    #q <- qr(fin.mat)
    #chas <- q$pivot[seq(q$rank)]
    #fin.mat <- as.matrix(fin.mat[,chas])
    #############
    return(fin.mat)
  }
  
  
  ## remove all markers or columns that are all missing data
  all.na <- apply(data,2,function(x){length(which(is.na(x)))/length(x)})
  bad.na <- which(all.na==1)
  if(length(bad.na) > 0){
    data <- data[,-bad.na]
  }
  #############################
  if(by.allele){ ####&&&&&&&&&&&&&&&&&&&&&& use zero.one function
    user.code <- apply(data[,c(1:(round(dim(data)[2]/20)))], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
    AA <- sum(user.code, na.rm = TRUE)/length(user.code)
    if(AA > .9){ # means user is using single letter
      rnd <- rownames(data)
      data <- apply(data,2,gbs.to.bisnp);#W2[1:5,1:5]
      rownames(data) <- rnd
    }
    M <- zero.one(data)
    
  }else{ ###&&&&&&&&&&&&&&&&&&&&&&&&
    n.g <- apply(data,2,function(x){length(table(x))})
    bad <- which(n.g > 3)
    if(length(bad) == dim(data)[2]){
      cat("Error. All your markers are multiallelic. This function requires at least one bi-allelic marker\n")
    }
    
    # tells you which markers have double letter code, i.e. TT instead of T
    # 1: has only one letter
    # 0: has two letters
    user.code <- apply(data[,c(1:(round(dim(data)[2]/20)))], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
    AA <- sum(user.code, na.rm = TRUE)/length(user.code)
    if(AA > .9){
      rrn <- rownames(data)
      
      cat("Converting GBS or single-letter code to biallelic code\n")
      if(silent){
        data <- apply(data, 2,gbs.to.bisnp)
      }else{
        data <- apply_pb(data, 2,gbs.to.bisnp) 
      }
      rownames(data) <- rrn
      data <- as.data.frame(data)
    }
    #### apply with progress bar ######
    s1 <- rownames(data)
    s2 <- colnames(data)
    data <- as.data.frame(t(data))
    rownames(data) <- s2
    colnames(data) <- s1
    bases <- c("A", "C", "G", "T","l","m","n","p","h","k","-","+","e","f","g","a","b","c","d")
    ## get reference allele function
    get.ref <- function(x, format) {
      if (format == "numeric") {
        ref.alt <- c(0, 1)
      }
      if (format == "AB") {
        ref.alt <- c("A", "B")
      }
      if (format == "ATCG") {
        y <- paste(na.omit(x), collapse = "")
        ans <- apply(array(bases), 1, function(z, y) {
          length(grep(z, y, fixed = T))
        }, y)
        if (sum(ans) > 2) {
          ref.alt <- (bases[which(ans == 1)])[1:2]
          #stop("Error in genotype matrix: More than 2 alleles")
        }
        if (sum(ans) == 2) {
          ref.alt <- bases[which(ans == 1)]
        }
        if (sum(ans) == 1) {
          ref.alt <- c(bases[which(ans == 1)], NA)
        }
      }
      return(ref.alt)
    }
    
    get.multi <- function(x, format) {
      if (format == "numeric") {
        ref.alt <- c(0, 1)
      }
      if (format == "AB") {
        ref.alt <- c("A", "B")
      }
      if (format == "ATCG") {
        y <- paste(na.omit(x), collapse = "")
        ans <- apply(array(bases), 1, function(z, y) {
          length(grep(z, y, fixed = T))
        }, y)
        if (sum(ans) > 2) {
          ref.alt <- TRUE
        }
        if (sum(ans) == 2) {
          ref.alt <- FALSE
        }
        if (sum(ans) == 1) {
          ref.alt <- FALSE
        }
      }
      return(ref.alt)
    }
    
    ####################################
    ## convert to matrix format
    ####################################
    markers <- as.matrix(data)
    ####################################
    # get reference alleles
    ####################################
    cat("Obtaining reference alleles\n")
    if(silent){
      tmp <- apply(markers, 1, get.ref, format=format)
    }else{
      tmp <- apply_pb(markers, 1, get.ref, format=format) 
    }
    
    if(multi){ # if markers with multiple alleles should be removed
      cat("Checking for markers with more than 2 alleles. If found will be removed.\n")
      if(silent){
        tmpo <- apply(markers, 1, get.multi, format = format)
      }else{
        tmpo <- apply_pb(markers, 1, get.multi, format = format) 
      }
      ###&&&&&&&&&&&& HERE WE MUST INSERT THE NEW FUNCTIONALITY, WHERE WE DETECTED MULTIPLE ALLELES
      multi.allelic <- which(!tmpo) # good markers
      markers <- markers[multi.allelic,]
      tmp <- tmp[, multi.allelic]
    }
    
    Ref <- tmp[1, ]
    Alt <- tmp[2, ]
    ####################################
    ## bind reference allele and markers and convert to numeric format based on the 
    # reference/alternate allele found
    ####################################
    cat("Converting to numeric format\n")
    if(silent){
      M <- apply(cbind(Ref, markers), 1, function(x) {
        y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
        ans <- as.integer(lapply(y, function(z) {
          ifelse(z[1] < 0, ploidy, ploidy - length(z))
        }))
        return(ans)
      })
    }else{
      M <- apply_pb(cbind(Ref, markers), 1, function(x) {
        y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
        ans <- as.integer(lapply(y, function(z) {
          ifelse(z[1] < 0, ploidy, ploidy - length(z))
        }))
        return(ans)
      })
    }
    
    gid.geno <- s1 #colnames(geno)
    rownames(M) <- gid.geno
    ####################################
    # identify bad markers
    ####################################
    bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
    if (bad > 0) {
      stop("Invalid marker calls.")
    }
    ####################################
    # by column or markers calculate MAF
    ####################################
    cat("Calculating minor allele frequency (MAF)\n")
    if(silent){
      MAF <- apply(M, 2, function(x) {
        AF <- mean(x, na.rm = T)/ploidy
        MAF <- ifelse(AF > 0.5, 1 - AF, AF)
      })
    }else{
      MAF <- apply_pb(M, 2, function(x) {
        AF <- mean(x, na.rm = T)/ploidy
        MAF <- ifelse(AF > 0.5, 1 - AF, AF)
      })
    }
    ####################################
    # which markers have MAF > 0, JUST GET THOSE
    ####################################
    polymorphic <- which(MAF > maf)
    M <- M[, polymorphic]
    ####################################
    # function to impute markers with the mode
    ####################################
    impute.mode <- function(x) {
      ix <- which(is.na(x))
      if (length(ix) > 0) {
        x[ix] <- as.integer(names(which.max(table(x))))
      }
      return(x)
    }
    # time to impute
    missing <- which(is.na(M))
    if (length(missing) > 0) {
      cat("Imputing missing data with mode \n")
      if(silent){
        M <- apply(M, 2, impute.mode)
      }else{
        M <- apply_pb(M, 2, impute.mode)
      }
    }
    if(ploidy == 2){
      M <- M - 1
    }
  }
  #rownames(M) <- rownames(data)
  ####################################
  return(M)
}