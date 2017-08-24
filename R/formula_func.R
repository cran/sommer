 at <- function(x, levs){ # how you handle x
#   # if(missing(levs)){
#   #   env = parent.frame()
#   #   levs <- unique(env$data[,deparse(x)])
#   # }
   dd <- model.matrix(~x - 1,data.frame(x))
   colnames(dd) <- substring(colnames(dd),2)
#   if(!missing(levs)){
     dd <- dd[,levs]
#   }
#   #colnames(dd) <- levs
#   return(dd)
 }

#at <- diag

g <- function(x){x}

#at <- function(x){x}

and <- function(x){x}

us <- function(x){x}

eig <- function(x){x}

grp <- function(x){x}

