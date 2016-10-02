manhattan <- function(map, col=NULL, fdr.level=0.05, show.fdr=TRUE, PVCN=NULL){
  
  if(!is.null(PVCN)){
    colnames(map)[which(colnames(map)==PVCN)] <- "p.val"
  }
  
  required.names <- c("Chrom","Position","p.val")
  
  che <- which(names(map)%in%required.names)
  if(length(che) < 3){
    stop("Column names; 'Chrom','Position' and 'p.val' need 
         to be present in the data frame provided.",call. = FALSE)
  }
  
  map <- map[with(map, order(Chrom, Position)), ]
  
  
  yylim <- ceiling(max(map$p.val))
  if(is.null(col)){
    col.scheme <- rep((transp(c("cadetblue","red"))),30)
  }else{
    col.scheme <- rep(col,30)
  }
  ffr <- fdr(map$p.val, fdr.level=fdr.level)$fdr.10
  plot(map$p.val, bty="n", col=col.scheme[factor(map$Chrom, levels = unique(map$Chrom, na.rm=TRUE))], xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), pch=20, cex=2, las=2, ylim=c(0,yylim))
  init.mrks <- apply(data.frame(unique(map$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=map$Chrom)
  fin.mrks <- apply(data.frame(unique(map$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=map$Chrom)
  inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
  axis(side=1, at=inter.mrks, labels=paste("Chr",unique(map$Chrom),sep=""), cex.axis=.5)
  if(show.fdr){
  abline(h=ffr, col="slateblue4", lty=3, lwd=2)
  legend("topright", legend=paste("FDR(",fdr.level,")=",round(ffr,2), sep=""), 
         bty="n", lty=3, lwd=2, col="slateblue4", cex=0.8)
  }
  #abline(v=marker, lty=3, col="red")
  #legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
  #marker <- as.character(map$Locus[marker])
}