
name.change <- function(x, separ="P", maxn=999){
  
  if(is.numeric(x)){
    stop("This function only works for character vectors.",call. = FALSE)
  }
  testo <- (gsub(separ,"",as.character(x)))
  
  testo2 <- as.numeric(testo)
  if(length(which(is.na(testo2)))==length(testo2)){
    #?name.change
    stop("When removing the separ argument in the vector doesn't yield a numeric vector. \nPlease check your vector 'x' and the 'separ' argument.",call. = FALSE)
  }
  
  x <- as.character(x)
  nosep <- gsub(separ,"",x)
  nosepn <- as.numeric(nosep)
  # posible divisions, which is closer to 1
  bb<- c(10,100,1000,10000,100000)
  base <- abs(10- maxn/bb)
  whereweare <- which(base == min(base))
  
  toadd <- gsub("1","",as.character(bb[1:whereweare]))
  toadd <- rev(toadd)
  ## get which jumps need to be considered and go from smaller to higher
  geto <- bb[1:whereweare]
  
  xn <- x
  
  for(u in 1:length(geto)){
    if(u ==1){
      l10 <- which(nosepn < 10)
      xn[l10] <- paste(separ,toadd[u],nosepn[l10],sep="")
    }else{
      l100 <- which(nosepn < geto[u] & nosepn >= geto[u-1])
      xn[l100] <- paste(separ,toadd[u],nosepn[l100],sep="")
    }
  }
  return(xn)
  
}