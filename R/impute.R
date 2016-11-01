imputev <- function(x, method="median"){
  if(is.numeric(x)){
    if(method=="mean"){
      x[which(is.na(x))] <- mean(x,na.rm=TRUE)
    }else if(method=="median"){
      x[which(is.na(x))] <- median(x,na.rm=TRUE)
    }else{
      x[which(is.na(x))] <- mean(x,na.rm=TRUE)
    }
  }else{
    if(method=="mean"){
      stop("Method 'mean' is not available for non-numeric vectors.",call. = FALSE)
    }else if(method=="median"){
      tt <- table(x)
      x[which(is.na(x))] <-  names(tt)[which(tt==max(tt))]
    }else{
      x[which(is.na(x))] <-  names(tt)[which(tt==max(tt))]
    }
  }
  return(x)
}