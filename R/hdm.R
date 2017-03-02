hdm <- function(dat){
  ss1 <- intersect(colnames(dat), c("female","male"))
  if(length(ss1)==0){
    stop("Columns 'female' and 'male' are rquired for building a hald diallel matrix",call.=FALSE)
  }
  dat2 <- dat[,ss1]; head(dat2)
  dat2$female <- as.factor(dat2$female)
  dat2$male <- as.factor(dat2$male)
  S1 <- model.matrix(~female-1, dat2)
  S2 <- model.matrix(~male-1, dat2)
  colnames(S1) <- gsub("female","",colnames(S1))
  colnames(S2) <- gsub("male","", colnames(S2))
  # make a matrix with all possible names full of zeros
  levo <- sort(unique(c(colnames(S1),colnames(S2)))); levo
  S3 <- matrix(0,nrow=dim(dat2)[1],ncol=length(levo))
  colnames(S3) <- levo
  S3[,colnames(S1)] <- S1 # ADD FEMALES
  S3[,colnames(S2)] <- S3[,colnames(S2)] + S2[,colnames(S2)] ## add males
  return(S3)
}
# 
# hdm <- function(data){
#   v1 <- data$female
#   v2 <- data$male
#   factoro <- FALSE
#   # condition to make sure is numeric
#   if((!is.numeric(v1) | !is.numeric(v2))){
#     #stop()
#     #cat("Please provide the data argument in numeric format")
#     factoro <- TRUE
#     ## names of matrix
#     vnames <- levels(as.factor(c(as.character(v1),as.character(v2))))
#     v3 <- as.numeric(as.factor(c(as.character(v1),as.character(v2))))
#     v4 <- as.factor(c(as.character(v1),as.character(v2)))
#     v1.1 <- v3[1:length(v1)]
#     v2.1 <- v3[(length(v1)+1):length(v3)]
#     v1 <- v1.1
#     v2 <- v2.1
#   }
#   
#   # vectors need to be numeric to work
#   ncol <- max(c(v1,v2), na.rm=TRUE)
#   nrow <- length(v1)
#   nana <- sort(unique(c(v1,v2)), decreasing=FALSE)
#   if(length(v1) != length(v2)){
#     stop
#     print("your maternal and paternal vectors need to have the same size")
#   }
#   Z <- matrix(0,nrow=nrow, ncol=ncol)
#   colnames(Z) <- nana
#   ## for parent 1 accomodate the 1's
#   for(i in 1:length(v1)){
#     ww1 <- which(colnames(Z) %in% v1[i])
#     Z[i,ww1] <- 1
#   }
#   ## for parent 2 accomodate the 1's
#   for(j in 1:length(v2)){
#     ww2 <- which(colnames(Z) %in% v2[j])
#     Z[j,ww2] <- 1
#   }
#   # if user provided factor instead of numeric
#   if(factoro){
#   colnames(Z) <- vnames
#   }
#   attributes(Z)$hdm <- TRUE
#   return(Z)
# }