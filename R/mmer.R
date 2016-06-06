mmer <- function(Y, X=NULL, Z=NULL, W=NULL, R=NULL, method="AI", REML=TRUE, iters=30, draw=FALSE, init=NULL, n.PC=0, P3D=TRUE, models="additive", ploidy=2, min.MAF=0.05, silent=FALSE, family=NULL, constraint=TRUE, sherman=FALSE, EIGEND=FALSE, Fishers=FALSE, gss=TRUE, forced=NULL, full.rank=TRUE, map=NULL, fdr.level=0.05, manh.col=NULL, gwas.plots=TRUE, n.cores=1,lmerHELP=FALSE){
  diso <- dim(as.data.frame(Y))[2]
  
  #########*****************************
  #########*****************************
  if(diso > 1){ # IF MULTIPLE RESPONSES
    NC <- detectCores() #parallel package
    
    #######$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #######$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if(n.cores > NC){ # USER USE MORE CORES THAN AVAILABLE
      stop(paste("You are selecting more cores than the",NC,"available in your computer"),call.=FALSE)
    }else{ # USER SELECTS LESS OR ENOUGH CORES THAN AVAILABLE
      
      ######^^^^^^^^^^^^^^^^^^^^^^^
      ######^^^^^^^^^^^^^^^^^^^^^^^
      if(NC >= n.cores){ # there's still more cores available than selected or enough
        
        if(diso > n.cores){ # and user has more responses than cores used
          cat(paste("Cores used in your computer:",n.cores, "of", NC,"\nFeel free to modify the 'n.core' argument to run each response in parallel.\nArguments 'draw', 'gwas.plots' will be set to FALSE")) 
          
          i=NULL
          cl <- makeCluster(n.cores) # Make clusters registerDoSNOW(cl)
          RES <- foreach (i=1:diso) %dopar% mmerSNOW(y=Y[,i], X=X, Z=Z, R=R, W=W, method=method, REML=REML, iters=iters, draw=FALSE, init=init, silent=TRUE, constraint=constraint, sherman=sherman, EIGEND=EIGEND, gss=gss, forced=forced, map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=FALSE,lmerHELP=lmerHELP)
          stopCluster(cl)
          names(RES) <- colnames(Y)
          class(RES)<-c("mmerM")
          
        }else{ # user has the same or less responses that the #of cores he is using
          
          cat(paste("Cores used in your computer:",n.cores, "for", diso,"responses.\nArguments 'draw', 'gwas.plots' will be set to FALSE"))
          i=NULL
          cl <- makeCluster(n.cores) # Make clusters registerDoSNOW(cl)
          RES <- foreach (i=1:diso) %dopar% mmerSNOW(y=Y[,i], X=X, Z=Z, R=R, W=W, method=method, REML=REML, iters=iters, draw=FALSE, init=init, silent=TRUE, constraint=constraint, sherman=sherman, EIGEND=EIGEND, gss=gss, forced=forced, map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=FALSE,lmerHELP=lmerHELP)
          stopCluster(cl)
          names(RES) <- colnames(Y)
          class(RES)<-c("mmerM")
          #RES$multi <- "M2"
        }
      }
    }
    #######$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #######$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    
  }else{ ## ONE RESPONSE
    RES <- mmerSNOW(y=Y, X=X, Z=Z, R=R, W=W, method=method, REML=REML, iters=iters, draw=draw, init=init, silent=silent, constraint=constraint, sherman=sherman, EIGEND=EIGEND, gss=gss, forced=forced, map=map, fdr.level=fdr.level, manh.col=manh.col,gwas.plots=gwas.plots,lmerHELP=lmerHELP)
    #RES$multi <- "M1"
  }
  #########*****************************
  #########*****************************
  
  ### "M1" is univariate
  ### "M2" is univariate run in parallel
  ### "M3" is multivariate run
  return(RES)
}


#### =========== ####
## SUMMARY FUNCTION #
#### =========== ####
"summary.mmer" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)
  groupss <- unlist(lapply(object$u.hat, function(y){dim(y)[1]}))
  groupss <- paste(groupss,collapse = " ")
  nn <- length(unlist(object$residuals))
  #cat("Linear mixed model fit by restricted maximum likelihood\n")
  #cat("=======================================================")
  #cat("\nMethod:")
  #print(object$method)
  method=object$method
  #cat("\n")
  
  logo <- c(logLik = object$LL, AIC=object$AIC, BIC=object$BIC)
  names(logo) = c("logLik","AIC","BIC")
  #print(c(logo), digits = digits)
  
  #cat("=======================================================")
  #cat("\nRandom effects:\n")
  w <- data.frame(VarianceComp = matrix(c(object$var.comp), ncol = 1))
  
  if(object$method == "EM"){
    row.names(w) <- rownames(object$var.comp)
  }
  if(object$method == "AI" | object$method == "EMMA" | object$method == "NR"){
    row.names(w) <- rownames(object$var.comp)
  }
  #print(w, digits = digits)
  #cat(paste("Number of obs:",nn," Groups:",groupss,"\n"))
  
  #cat("=======================================================")
  #cat("\nFixed effects:\n")
  coef <- data.frame(Value = as.matrix(object$beta.hat), Std.Error=(matrix(sqrt(diag(object$Var.beta.hat)),ncol=1)), t.value=(matrix((object$beta.hat-0)/sqrt(diag(object$Var.beta.hat)), ncol=1)))
  if(dim(coef)[1] == 1){rownames(coef) <- "Intercept"}
  #printCoefmat((coef))
  #cat("=======================================================")
  #cat("\nVar-Cov for Fixed effects (diagonals are variances):\n")
  varbhat <- data.frame(as.matrix(object$Var.beta.hat))
  rownames(varbhat) <- colnames(object$X)
  colnames(varbhat) <- rownames(varbhat) 
  if(dim(varbhat)[1] == 1){rownames(varbhat) <- "Intercept"}
  #printCoefmat(varbhat)
  #cat("=======================================================")
  #cat("\nInformation contained in this fitted model: \n* Variance components\n* Residuals and conditional residuals\n* Inverse phenotypic variance(V)\* BLUEs and BLUPs\n* Variance-covariance matrix for fixed effects\n* Variance-covariance matrix for random effects\n* Predicted error variance (PEV)\n* LogLikelihood\n* AIC and BIC\n* Fitted values\nUse the 'str' function to access such information")
  output <- list(groupss=groupss, nn=nn, logo=logo, w=w, coef=coef, varbhat=varbhat, method=method)
  attr(output, "class")<-c("summary.mmer", "list")
  return(output)
}

"print.summary.mmer"<-function (x, digits = max(3, getOption("digits") - 3),  ...){
  digits = max(3, getOption("digits") - 3)
  #groupss <- unlist(lapply(x$u.hat, function(y){dim(y)[1]}))
  #groupss <- paste(groupss,collapse = " ")
  #nn <- length(unlist(x$residuals))
  cat("\nInformation contained in this fitted model: \n* Variance components\n* Residuals and conditional residuals\n* BLUEs and BLUPs\n* Inverse phenotypic variance(V)\n* Variance-covariance matrix for fixed effects\n* Variance-covariance matrix for random effects\n* Predicted error variance (PEV)\n* LogLikelihood\n* AIC and BIC\n* Fitted values\nUse the 'str' function to access such information\n")
  cat("\n=======================================================")
  cat("\nLinear mixed model fit by restricted maximum likelihood\n")
  cat("********************  sommer 1.8  *********************\n")
  cat("=======================================================")
  cat("\nMethod:")
  print(x$method)
  cat("\n")
  
  #logo <- c(logLik = x$LL, AIC=x$AIC, BIC=x$BIC)
  #names(logo) = c("logLik","AIC","BIC")
  print(c(x$logo), digits = digits)
  
  cat("=======================================================")
  cat("\nRandom effects:\n")
  #w <- data.frame(VarianceComp = matrix(c(x$var.comp), ncol = 1))
  
  #if(x$method == "EM"){
  #  row.names(w) <- names(x$var.comp)
  #}
  #if(x$method == "AI" | x$method == "EMMA"){
  #  row.names(w) <- rownames(x$var.comp)
  #}
  print(x$w, digits = digits)
  cat(paste("Number of obs:",x$nn," Groups:",x$groupss,"\n"))
  
  cat("=======================================================")
  cat("\nFixed effects:\n")
  #coef <- data.frame(Value = x$beta.hat, Std.Error=(matrix(sqrt(diag(x$Var.beta.hat)),ncol=1)), t.value=(matrix((x$beta.hat-0)/sqrt(diag(x$Var.beta.hat)), ncol=1)))
  #if(dim(coef)[1] == 1){rownames(coef) <- "Intercept"}
  printCoefmat((x$coef))
  cat("=======================================================")
  cat("\nVar-Cov for Fixed effects:\n(diagonals are variances)\n")
  #varbhat <- data.frame(x$Var.beta.hat)
  #colnames(varbhat) <- rownames(varbhat) 
  #if(dim(varbhat)[1] == 1){rownames(varbhat) <- "Intercept"}
  printCoefmat(x$varbhat)
  cat("=======================================================")
  cat("\nUse the 'str' function to access all information\n\n")
}
#### =========== ####
## SUMMARY FUNCTION 222222222222223 #
#### =========== ####
"summary.mmerM" <- function(object, ...) {
  
  
    digits = max(3, getOption("digits") - 3)
    #forget <- length(object)
    
    groupss.nn <- cbind(
      do.call("rbind",lapply(object,function(x){lapply(x$u.hat, function(y){dim(y)[1]})})),
      do.call("rbind",lapply(object,function(x){length(unlist(x$residuals))}))
    ); 
    groupss.nn <- t(groupss.nn);
    
    if(is.null(rownames(groupss.nn))){
      re <- 1:(dim(groupss.nn)[1] - 1) # random effects
      rownames(groupss.nn) <- c(paste("Groups.u",re,sep="."),"No.obs")
    }else{
      rownames(groupss.nn)[dim(groupss.nn)[1]] <- c("No.obs")
    }
    
    
    
    method=object[[1]]$method
    #cat("\n")
    
    LLAIC <- rbind(do.call("cbind",lapply(object,function(x){x$LL})),
    do.call("cbind",lapply(object,function(x){x$AIC})),
    do.call("cbind",lapply(object,function(x){x$BIC}))
    )
    
    rownames(LLAIC) = c("logLik","AIC","BIC")
    
    w <- do.call("cbind",lapply(object,function(x){x$var.comp}))
    colnames(w) <- names(object)
    output <- list(groupss=groupss.nn, logo=LLAIC, w=w, method=method)
    attr(output, "class")<-c("summary.mmerM", "list")
  return(output)
}

"print.summary.mmerM"<-function (x, digits = max(3, getOption("digits") - 3),  ...){
  
    digits = max(3, getOption("digits") - 3)
    cat("Information contained in this structure: \n* Individual results for each response model\nDisplayed: \n* AIC and BIC summaries\n* Variance component summaries\nUse the '$' sign to access individual models\n")
    cat("=======================================================")
    cat("\nLinear mixed model fit by restricted maximum likelihood\n")
    cat("********************  sommer 1.8  *********************\n")
    cat("=======================================================")
    cat("\nMethod:")
    print(x$method)
    cat("")
    print(x$logo)#, digits = digits)
    cat("=======================================================")
    cat("\nVariance components:\n")
    print(x$w, digits = digits)
    cat("=======================================================")
    cat("\nGroups and observations\n")
    print(x$groupss, digits = digits)
    cat("=======================================================")
    cat("\nUse the '$' sign to access individual models\nYou can use 'summary' on individual models\nArguments set to FALSE for multiple models:\n'draw' and 'gwas.plots'\n")
}

#### =========== ######
## RESIDUALS FUNCTION #
#### =========== ######
"residuals.mmer" <- function(object, type="conditional", ...) {
  digits = max(3, getOption("digits") - 3)
  if(type=="conditional"){
    output<-round(object$cond.residuals,digits)
  }else{
    output<-round(object$residuals,digits)
  }
  return(output)
}

"print.residuals.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}


"residuals.mmerM" <- function(object, type="conditional", ...) {
  digits = max(3, getOption("digits") - 3)

  if(type=="conditional"){
    output<-do.call("cbind",lapply(object,function(x){x$cond.residuals}))
    colnames(output) <- names(object)
  }else{
    output<-do.call("cbind",lapply(object,function(x){x$residuals}))
    colnames(output) <- names(object)
  }
  return(output)
}

"print.residuals.mmerM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

#### =========== ######
## RANEF FUNCTION #
#### =========== ######

"randef" <- function(object) {
  if(class(object)=="mmerM"){
    digits = max(3, getOption("digits") - 3)
    output<-do.call("cbind",lapply(object,function(x){do.call("rbind",x$u.hat)}))
    colnames(output) <- names(object)
  }else if(class(object)=="mmer"){
    digits = max(3, getOption("digits") - 3)
    output <- object$u.hat
  }

  return(output)
}

#"print.ranef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
#  print((x))
#}

#### =========== ######
## FIXEF FUNCTION #
#### =========== ######
"fixef.mmer" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)
  output <- object$beta.hat
  return(output)
}

"print.fixef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

"fixef.mmerM" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)
  output<-do.call("cbind",lapply(object,function(x){x$beta.hat}))
  colnames(output) <- names(object)
  return(output)
}

"print.fixef.mmerM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

#### =========== ####
## FITTED FUNCTION ##
#### =========== ####
"fitted.mmer" <- function(object, type="complete", ...) {
  #type="complete" 
  digits = max(3, getOption("digits") - 3)
  if(type=="complete"){
    round(object$fitted.y,digits)
  }else{
    round(object$fitted.u,digits)
  }
}

"print.fitted.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 

"fitted.mmerM" <- function(object, type="complete", ...) {
  #type="complete" 
  digits = max(3, getOption("digits") - 3)
  if(type=="complete"){
    output<-do.call("cbind",lapply(object,function(x){x$fitted.y}))
    colnames(output) <- names(object)
  }else{
    output<-do.call("cbind",lapply(object,function(x){x$fitted.u}))
    colnames(output) <- names(object)
  }
  return(output)
}

"print.fitted.mmerM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 
#### =========== ####
## COEF FUNCTION ####
#### =========== ####
"coef.mmer" <- function(object, ...){
  object$beta.hat
}

"print.coef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 

"coef.mmerM" <- function(object, ...){
  output<-do.call("cbind",lapply(object,function(x){x$beta.hat}))
  colnames(output) <- names(object)
  nana <- colnames(object[[1]]$X)
  if(is.null(nana) & (dim(object[[1]]$X)[2] == 1)){
    rownames(output) <- "Intercept"
  }else{
    rownames(output) <- nana
  }
  return(output)
}

"print.coef.mmerM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 

#### =========== ####
## ANOVA FUNCTION ###
#### =========== ####
anova.mmer <- function(object, object2=NULL, ...) {
  signifo <- function(x){
    if(x >= 0 & x < 0.001){y="***"}
    if(x >= 0.001 & x < 0.01){y="**"}
    if(x >= 0.01 & x < 0.05){y="*"}
    if(x >= 0.05 & x < 0.1){y="."}
    if(x > 0.1){y=""}
    return(y)
  }
  ########################################
  digits = max(3, getOption("digits") - 3)
  if(is.null(object2)){
    stop("The 'anova' function for the sommer package only works to compare mixed models by likelihood ratio tests (LRT), was not intended to provide regular sum of squares output.")
  }else{
    #if(object$maxim){ # user used REML=TRUE, not possible to do LRT
    #  stop("Please fit the models using ML instead of REML by setting the argument REML=FALSE and try again")
    #}else{ #if user used REML=FALSE, then proceed
    if(object$method != object2$method){
      stop("Error! When comparing models please use the same method for the fitted models.")
    }else{
      dis=c(dim(as.matrix(object$var.comp))[1]+dim(object$beta.hat)[1],
            dim(as.matrix(object2$var.comp))[1]+dim(object2$beta.hat)[1]) # dimensions
      mods=c("mod1","mod2")
      lls=c(object$LL, object2$LL) # likelihoods
      aics=c(object$AIC, object2$AIC) # AIC's
      bics=c(object$BIC, object2$BIC) # AIC's
      vv=which(dis == max(dis))[1] # which has more variance components BIGGER
      vv2=c(1:2)[which(c(1:2)!= vv)] # SMALLER
      LR = (lls[vv] - lls[vv2])
      r.stat= abs(-2*((LR))) # -2(LL1 - LL2)
      df=dis[vv]-dis[vv2]
      chichi=round(pchisq((r.stat), df, lower.tail=FALSE),16)
      chichi2=paste(as.character((chichi)),signifo(chichi), sep=" ")
      ### construct the table
      cat("Likelihood ratio test for mixed models\n")
      cat("==============================================================\n")
      result=data.frame(Df=c(dis[vv],dis[vv2]), AIC=c(aics[vv],aics[vv2]), 
                        BIC=c(bics[vv],bics[vv2]), loLik=c(lls[vv],lls[vv2]), 
                        Chisq=c("",as.character(round(r.stat,5))), 
                        ChiDf=c("",as.character(df)), PrChisq=c("",as.character(chichi2 )))
      rownames(result) <- c(mods[vv],mods[vv2])
      print(result)
      cat("==============================================================\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    }
    #}
  }
}

anova.mmerM <- function(object, object2=NULL, ...) {
  cat("'anova' function only works for individual models to perform likelihood ratio tests (LRT).\nYou can access individuals models using the '$' sign and then use 'anova'.")
}
#### =========== ####
## PLOTING FUNCTION #
#### =========== ####
plot.mmer <- function(x, ...) {
  digits = max(3, getOption("digits") - 3)
  transp <- function (col, alpha = 0.5){
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,c[3]/255, alpha))
    return(res)
  }
  # std vs residuals, QQplot (std vs teor quantiles), sqrt(std residuals) vs fitted, std res vs leverage = cook's distance
  layout(matrix(1:4,2,2))
  plot(x$fitted.y.good, scale(x$cond.residuals), pch=20, col=transp("cadetblue"), ylab="Std Residuals", xlab="Fitted values", main="Residual vs Fitted", bty="n", ...); grid()
  plot(x$fitted.y.good, sqrt(abs((scale(x$cond.residuals)))), pch=20, col=transp("thistle4"), ylab="Sqrt Abs Std Residuals", xlab="Fitted values", main="Scale-Location",bty="n", ...);grid()
  qqnorm(scale(x$cond.residuals), pch=20, col=transp("tomato1"), ylab="Std Residuals", bty="n",...); grid()
  hat <- x$X%*%solve(t(x$X)%*%x$V.inv%*%x$X)%*%t(x$X)%*%x$V.inv # leverage including variance from random effects H= X(X'V-X)X'V-
  plot(diag(hat), scale(x$cond.residuals), pch=20, col=transp("springgreen3"), ylab="Std Residuals", xlab="Leverage", main="Residual vs Leverage", bty="n", ...); grid()
  #####################
  layout(matrix(1,1,1))
}

plot.mmerM <- function(x, ...) {
  cat("'plot' function only works for individual models to check diagnostic plots.\nYou can access individuals models using the '$' sign and then use 'plot'.")
}
##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "2.1"))
    stop("This package requires R 2.1 or later")
  assign(".sommer.home", file.path(library, pkg),
         pos=match("package:sommer", search()))
  sommer.version = "1.8 (2016-06-06)"
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### check which version is more recent
  #yyy <- 1.8
  #chooseCRANmirror(ind=114)
  #xxx <- available.packages(contriburl = contrib.url(repos="http://mirror.las.iastate.edu/CRAN/", type = getOption("pkgType")))
  
  #xxx <- available.packages()
  #current <- as.numeric(xxx["sommer","Version"])
  ### final check
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  assign(".sommer.version", sommer.version, pos=match("package:sommer", search()))
  if(interactive())
  {
    packageStartupMessage(paste("## ========================================================= ## "),appendLF=TRUE)
    packageStartupMessage(paste("# Solving Mixed Model Equations in R (sommer) ", sommer.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage(paste("# Mixed models allowing covariance structures in random effects"),appendLF=TRUE)
    packageStartupMessage("# Author: Giovanny Covarrubias-Pazaran",appendLF=TRUE)
    packageStartupMessage("# Published: PLOS ONE",appendLF=TRUE)
    packageStartupMessage("# Supported by the Council of Science and Technology (CONACYT)", appendLF=TRUE)
    packageStartupMessage("# Type 'help(sommer)' for summary information",appendLF=TRUE)
    packageStartupMessage("# Type 'vignette('sommer')' for a short tutorial",appendLF=TRUE)
    packageStartupMessage("# Type 'citation('sommer')' to know how to cite sommer",appendLF=TRUE)
    packageStartupMessage(paste("## ========================================================= ## "),appendLF=TRUE)
    packageStartupMessage("PLEASE UPDATE 'sommer' EVERY MONTH FOR NEW FEATURES USING",appendLF=TRUE)
    packageStartupMessage(paste("install.packages('sommer') -- your version: ", sommer.version, ". ",sep=""),appendLF=TRUE)
    
    #if(yyy > current){ # yyy < current in CRAN
    #  packageStartupMessage(paste("Version",current,"is now available."),appendLF=TRUE) # version current
    #  packageStartupMessage(paste("Please update 'sommer' installing the new version."),appendLF=TRUE) # version current
    #}
    
  }
  invisible()
}

#.onLoad <- function(library, pkg) {
#  data(x)
#  library(audio)
#  packageStartupMessage(play(x))
#}