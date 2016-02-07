fdr <- function(p, fdr.level=0.05){

  pval <- 10^-p

  pvalA <- p.adjust(pval, method="fdr")

  
  pvalA.l10 <- -log(pvalA, 10)

  fdr.ad <- -log(fdr.level, 10)

  vv <- which(pvalA.l10 > fdr.ad)
  if(length(vv)>0){
    fdr.or <- min(p[vv])
  }else{
    fdr.or <- NULL
  }

  result <- list(p.ad=pvalA.l10, fdr=fdr.ad, p.or=p, fdr.or=fdr.or )
  return(result)
}