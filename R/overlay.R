# overlay <-function(dat, r1=1, r2=1){
#   dat <- as.data.frame(dat)
#   if(dim(dat)[2]>2){
#     stop("Please provide a dataframe with 2 columns to do the overlay\n", call. = FALSE)
#   }
#   ss1 <- colnames(dat)
#   dat2 <- dat[,ss1]; head(dat2)
#   
#   fem <- ss1[1]
#   mal <- ss1[2]
#   dat2[,fem] <- as.factor(dat2[,fem])
#   dat2[,mal]<- as.factor(dat2[,mal])
#   S1 <- model.matrix(as.formula(paste("~",fem,"-1")), dat2)
#   S2 <- model.matrix(as.formula(paste("~",mal,"-1")), dat2)
#   colnames(S1) <- gsub(fem,"",colnames(S1))
#   colnames(S2) <- gsub(mal,"", colnames(S2))
#   # make a matrix with all possible names full of zeros
#   levo <- sort(unique(c(colnames(S1),colnames(S2)))); levo
#   S3 <- matrix(0,nrow=dim(dat2)[1],ncol=length(levo)); rownames(S3) <- rownames(dat2)
#   colnames(S3) <- levo
#   S3[rownames(S1),colnames(S1)] <- S1*r1 # ADD FEMALES
#   S3[rownames(S2),colnames(S2)] <- S3[rownames(S2),colnames(S2)] + (S2[rownames(S2),colnames(S2)]*r2) ## add males
#   return(S3)
# }


overlay <-function(dat, rlist=NULL, prefix=NULL){
  
  if(is.null(dim(dat))){
    stop("Please provide a data frame to the overlay function, not a vector.\n", call. = FALSE)
  }
  dat <- as.data.frame(dat)
  if(is.null(rlist)){
    #rlist <- vector(mode="list",length = dim(dat)[2])
    rlist <- as.list(rep(1, dim(dat)[2]))
  }
  
  ss1 <- colnames(dat)
  dat2 <- as.data.frame(dat[,ss1]); head(dat2)
  colnames(dat2) <- ss1
  
  femlist <- list()
  S1list <- list()
  ## convert to factor such columns
  ## and store matrices
  for(i in 1:length(ss1)){
    femlist[[i]] <- ss1[i]
    dat2[,femlist[[i]]] <- as.factor(dat2[,femlist[[i]]])
    S1 <- model.matrix(as.formula(paste("~",femlist[[i]],"-1")), dat2)
    colnames(S1) <- gsub(femlist[[i]],"",colnames(S1))
    S1list[[i]] <- S1
  }
  # make a matrix with all possible names full of zeros
  levo <- sort(unique(unlist(lapply(S1list, function(x){colnames(x)}))))
  S3 <- matrix(0,nrow=dim(dat2)[1],ncol=length(levo)); rownames(S3) <- rownames(dat2)
  colnames(S3) <- levo
  for(i in 1:length(S1list)){
    if(i==1){
      S3[rownames(S1list[[i]]),colnames(S1list[[i]])] <- S1list[[i]]*rlist[[i]] # ADD FEMALES
    }else{
      S3[rownames(S1list[[i]]),colnames(S1list[[i]])] <- S3[rownames(S1list[[i]]),colnames(S1list[[i]])] + (S1list[[i]][rownames(S1list[[i]]),colnames(S1list[[i]])]*rlist[[i]]) ## add males
    }
  }
  if(!is.null(prefix)){
    colnames(S3) <- paste(prefix,colnames(S3),sep="")
  }
  return(S3)
}