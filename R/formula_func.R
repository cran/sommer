at <- function(x, levs){ # how you handle x
  dd <- model.matrix(~x - 1,data.frame(x))
  colnames(dd) <- substring(colnames(dd),2)
  dd <- dd[,levs]
  return(dd)
}

g <- function(x){x}


and <- function(x){x}

us <- function(x){x}