spatPlots <- function(object, by=NULL, colfunc=NULL,row="ROW",range="RANGE", wire=FALSE){
  
  dat <- object$data
  
  if(is.null(colfunc)){
    colfunc <- colorRampPalette(c("gold","springgreen","steelblue4"))
  }
  if(missing(dat)){
    stop("dat argument needs to be provided.\n",call. = FALSE)
  }
  if(missing(object)){
    stop("asreml object needs to be provided.\n",call. = FALSE)
  }
  rr <- which(colnames(dat) %in% c(row,range))
  if(length(rr) != 2){
    stop("row and range arguments didn't match column names in your model.\n", call. = FALSE)
  }
  
  #object$call
  
  refs <- apply(data.frame(object$random.effs),1,function(x){
    pp <- strsplit(x,":")[[1]]
    return(pp[length(pp)])
  })
  (refs2 <- unique(refs))
  ## find which index belong to each random effect
  index <- numeric()
  for(u in 1:length(refs2)){
    v <- which(refs == refs2[u])
    index[v] <- u
  }
  ## correct in case if there's 2 dimnesional splines join them
  v <- grep("_2Dspl",refs)
  if(length(v)>0){
    index[v] <- min(index[v])
    refs[v] <- "spl2D"
  }
  refs2 <- unique(refs)
  ##
  uindex <- unique(index)
  
  fits <- list()
  for(ui in uindex){ # ui <- uindex[1]
    we <- which(index==ui)
    fits[[ui]] <- apply(do.call(cbind,object$Zus[we]),1,sum)
  }
  fits[[length(fits)+1]] <- object$res.ordim
  fits <- as.data.frame(do.call(cbind,fits))
  refs2 <- gsub("\\(",".",refs2)
  refs2 <- gsub(")",".",refs2)
  colnames(fits) <- c(paste("fit",refs2,sep="_"),"residuals")
  
  fits <- data.frame(dat,fits)
  colos <- c(paste("fit",refs2,sep="_"),"residuals")
  
  for(u in 1:length(colos)){
    
    resp <- colos[u]
    
    if(is.null(by)){
      form <- as.formula(paste(resp,"~",row,"*",range))
    }else{
      form <- as.formula(paste(resp,"~",row,"*",range,"|",by))
    }
    
    if(wire){ # wireframe
      print(wireframe(form, data=fits,  
                      aspect=c(61/87,0.4), drape=TRUE,
                      light.source=c(10,0,10), #=colfunc,
                      main=resp)
      )
    }else{ # levelplot
      print(levelplot(form, data=fits, main=resp, col.regions = colfunc))
    }
    
  }
  #head(fittedv)
  return(fits)
}