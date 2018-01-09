variogram <- function (x, xcoor="ROW", ycoor="RANGE", zcoor=NULL, by=NULL, ...) {
  UseMethod("variogram")
}

variogram.MMERM <- function (x, xcoor="ROW", ycoor="RANGE", zcoor=NULL, by=NULL, ...){
  #x is a dataset with columns Residuals, xcoor, ycoor
  # xcoor is the name of column for the xcoordinate
  # ycoor is the name of column for the ycoordinate
  #library(data.table)
  
  x0 <- x$data
  
  if(is.null(by)){
    x0$FIELDINST <- "F1"
    x0$FIELDINST <- as.factor(x0$FIELDINST)
    by="FIELDINST"
  }else{
    v <- which(colnames(x0) == by)
    if(length(v)==0){stop("by argument not found in the column names of x0", call. = FALSE)}
    x0[,by] <- as.factor(x0[,by])
  }
  
  if(is.null(zcoor)){
    zcoor <- "Residuals"
    x0[,zcoor] <- x$res.ordim
  }else{
    www <- which(names(x$Zus) %in% zcoor)
    if(length(www)==1){
      x0[,zcoor] <- x$Zus[[zcoor]]
    }else if(length(www)>1){
      zcoor <- paste(zcoor, collapse = "_")
      gggg <- Reduce("+",x$Zus[www])
      x0[,zcoor] <- gggg
    }else{
      stop("zcoor not found in your model", call. = FALSE)
    }
  }
  
  are <- which(colnames(x0) %in% c(xcoor,ycoor,zcoor))
  if(length(are) < 3){
    stop("One or more of xcoor, ycoor, zcoor don't match your data frame.", call. = FALSE)
  }
  
  ## now add the zcoor
  
  x1<- split(x0, f=x0[,by])
  
  multires <- lapply(x1, function(x){
    x.coord <- x[, xcoor]
    y.coord <- x[, ycoor]
    residuals <- x[, zcoor]
    columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
    rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
    xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
    setkeyv(xy.coord, c("rows", "columns"))
    df <- data.table(columns = x.coord, rows = y.coord, residuals = residuals)
    setkeyv(df, c("rows", "columns"))
    df <- df[xy.coord]
    df <- df[order(df$columns, df$rows), ]
    resdiff <- c(outer(df$residuals, df$residuals, function(x, 
                                                            y) 0.5 * (x - y)^2))
    coldiff <- c(outer(df$columns, df$columns, function(x, y) abs(x - 
                                                                    y)))
    coldiff.u <- unique(coldiff)
    rowdiff <- c(outer(df$rows, df$rows, function(x, y) abs(x - 
                                                              y)))
    rowdiff.u <- unique(rowdiff)
    subsets <- split(resdiff, f = list(coldiff, rowdiff))
    value <- sapply(subsets, mean, na.rm = TRUE)
    length <- sapply(subsets, function(x) sum(!is.na(x)))
    length[-1] <- length[-1]/2
    res <- list(data = data.frame(value = value, length = length), 
                col.displacement = coldiff.u, row.displacement = rowdiff.u)
    class(res) <- "variogram.MMERM"
    return(res)
  })
  #print(str(multires))
  # multires <- as.data.frame(do.call(rbind,multires))
  # class(multires) <- "variogram.sommer"
  # 
  # funny <- as.formula(paste(zcoor,"~",xcoor,"*",ycoor,"|",by, sep="" ))
  # multires2 <- list(data=multires,funny=funny)
  class(multires) <- "variogram.MMERM"
  
  return(multires)
}

###########################
###########################
###########################
###########################
###########################

plot.variogram.MMERM <- function(x, stnd=TRUE, ...) {
  
  #for(u in 1:length(x)){
    
    x0 <- x#[[u]]
    min.length = 30
    if(stnd){x0$data$value <- scale(x0$data$value)}
    values <- matrix(replace(x0$data$value, x0$data$length < min.length,
                             NA), ncol = length(x0$col.displacement), nrow = length(x0$row.displacement),
                     byrow = TRUE)
    
    print(wireframe(values,drape=TRUE, #main=names(x)[u],
                    aspect = c(61/87, 0.4), #colorkey=TRUE,
                    light.source = c(10,0,10), #shade=TRUE,
                    #screen = list(x0 = -60, y = 50, z=20),
                    col.regions = topo.colors(100))
    )
    # light.source = c(0,0,10),
    # #region = TRUE,
    # col.regions = terrain.colors(100),
    # screen = list(x0 = -60, y = 50, z=20)))
    
    # persp(x0, y, z, theta = 135, phi = 30, col = colorRampPalette(c("blue", "pink"))(9500), scale = FALSE,
    #       ltheta = -120, shade = 0.75, border = NA, box = FALSE)
    
  #}
  # plot3Drgl::persp3Drgl(x0$row.displacement, x0$col.displacement,
  #                       values, xlab = "Row displacement", ylab = "Col displacement",
  #                       zlab = "", ticktype = "detailed", col = jet.colors(100))
}