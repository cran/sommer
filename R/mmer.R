mmer <- function(Y,X=NULL,Z=NULL,R=NULL,W=NULL,method="NR",init=NULL,iters=20,tolpar=1e-3,
                 tolparinv=1e-6,draw=FALSE,silent=FALSE, constraint=TRUE, 
                 EIGEND=FALSE, forced=NULL, IMP=FALSE, complete=TRUE, 
                 check.model=TRUE, restrained=NULL, REML=TRUE, init.equal=TRUE,
                 date.warning=TRUE){
  #R=NULL
  if (inherits(Y, "formula")){
    stop("\nYou have to use mmer2 function for formula-based models.\n", call. = FALSE)
  }
  diso <- dim(as.data.frame(Y))[2]
  ## control for 2-level list structure
  if(!is.list(Z)){
    stop("Please provide the Z parameter as a 2 level list structure.\nFor example for 2 random effects 'A' and 'B' do:\n    ETA <- list( A=list( Z=myZ1, K=myK1 ) , B=list( Z=myZ2, K=myK2 ) )\n    mod <- mmer(Y=y, Z=ETA)\nwhere Z's and K's are the incidence and var-covar matrices respectively.\nIf any Z or K is not provided, an identity matrix will be assumed. ",call. = FALSE)
  }else{
    if(!is.list(Z[[1]])){
      stop("Please provide the Z parameter as a 2 level list structure.\nFor example for 2 random effects 'A' and 'B' do:\n    ETA <- list( A=list( Z=myZ1, K=myK1 ) , B=list( Z=myZ2, K=myK2 ) )\n    mod <- mmer(Y=y, Z=ETA)\nwhere Z's and K's are the incidence and var-covar matrices respectively.\nIf any Z or K is not provided, an identity matrix will be assumed. ",call. = FALSE)
    }
  }
  ## control for Z-K names
  zzzkkk <- unlist(lapply(Z,function(x){length(names(x))}))
  badRE <- which(zzzkkk==0) # BAD RE WITH NO NAMES
  if(length(badRE)>0){
    stop("Please when specifying a random effect use the names; \n'Z' for incidence and 'K' for variance-covariance matrices.\nFor example for 1 random effect (i.e. named 'A') model do:\n    ETA <- list( A=list( Z=M1, K=M2) )\n    mod <- mmer(Y=y, Z=ETA)\nSpecifying at least one; Z or K. You need to specify if is a 'Z' or 'K' \nsince this is the only way the program distinguishes between matrices.",call. = FALSE)
  }
  
  my.year <- 2018
  my.month <- 6 #month when the user will start to get notifications the 1st day of next month
  ### if my month = 3, user will start to get notification in april 1st (next month)
  datee <- Sys.Date()
  year.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
  your.year <- year.mo.day[1]
  your.month <- year.mo.day[2]
  ## if your month is greater than my month you are outdated
  if(date.warning){
    if(your.month > my.month & your.year >= my.year){
      # error if your month is greater and your year is smaller
      cat("Version out of date. Please update sommer to the newest version using:\ninstall.packages('sommer') in a new session\n")
    }
  }
  #########*****************************
  ## make sure user don't provide the same names for random effects
  his.names <- names(Z)
  if(!is.null(his.names)){
    badnames <- which(duplicated(his.names))
    if(length(badnames)>0){
      his.names[badnames] <- paste(his.names[badnames],1:(length(his.names[badnames])), sep=".")
      names(Z) <- his.names
    }
  }
  dZ <- unlist(lapply(Z,function(x){dim(x$Z)[1]}))
  
  if(is.null(dZ)){ #sometimes user don't specify the Z matrices
    dZ <- unlist(lapply(Z,function(x){dim(x$K)[1]}))
  }
  if(!is.null(X)){
    dZ <- c(dZ,dim(X)[1])
  }
  dall <- unlist(c(dZ,dim(as.matrix(Y))[1]))
  if(length(which(!duplicated(dall))) > 1){
    if(is.null(X)){
      stop("Matrices Y and Z's should have the same number of individuals. \nPlease check the dimensions of your matrices.", call. = FALSE)
    }else{
      stop("Matrices Y, X and Z's should have the same number of individuals. \nPlease check the dimensions of your matrices.", call. = FALSE)
    }
  }
  #########*****************************
  for(bb in 1:length(Z)){
    ss1 <- colnames(Z[[bb]]$Z) == colnames(Z[[bb]]$K)
    if(length(which(!ss1))>0){
      print(paste("Names of columns in matrices Z and K for the",bb,"th random effect do not match.")) 
      print("This can lead to incorrect estimation of variance components. Double check.")
    }
  }
  #########*****************************
  #   if(!is.null(W)){
  #     cat("Response is imputed for estimation of variance components in GWAS models.\n")
  #   }
  #########*****************************
  
  if(!is.null(X)){
    if(is.list(X)){
      stop("Multivariate models only accept one incidence matrix for fixed effects (X). Please modifiy your X argument.",call. = FALSE)
    }
  }
  #if(!silent){cat("Running multivariate model\n")}
  
  if(check.model){
    if(is.list(Z)){
      if(is.list(Z[[1]])){ ### -- if is a 2 level list -- ##
        provided <- lapply(Z, names)
        for(s in 1:length(provided)){ #for each random effect =============================
          provided2 <- names(Z[[s]])
          if(length(provided2) ==1){ #----the 's' random effect has one matrix only----
            if(provided2 == "K"){ #user only provided K
              #zz <- diag(length(y))#model.matrix(~rownames(Z[[s]][[1]]))
              zz <- diag(nrow(as.matrix(Y)))
              #colnames(zz) <- rownames(Z[[s]][[1]])
              Z[[s]] <- list(Z=zz, K=Z[[s]][[1]])
            }
            if(provided2 == "Z"){ # user only provided Z
              #kk <- diag(dim(Z[[s]][[1]])[2])
              kk <- diag(dim(Z[[s]][[1]])[2])
              attributes(kk)$diagon <- TRUE
              #rownames(kk) <- colnames(Z[[s]][[1]]); colnames(kk) <- rownames(kk)
              Z[[s]] <- list(Z=Z[[s]][[1]],K=kk) 
            }
          }else{ #----the 's' random effect has two matrices----
            dido<-lapply(Z[[s]], dim) # dimensions of Z and K
            condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
            # condition, column size on Z matches with a square matrix K
            if(!condi){
              cat(paste("ERROR! In the",s,"th random effect you have provided or created an incidence \nmatrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2]," matrix \nas a variance-covariance matrix. Please double check your matrices.")
              stop()
            }
          }#---------------------------------------------------------------------------
        } #for each random effect end =================================================
      }else{ # if is a one-level list !!!!!!!!!!!!!
        if(length(Z) == 1){ ## -- if the user only provided one matrix -- ##
          provided <- names(Z)
          if(provided == "K"){
            #zz <- diag(length(y))
            zz <- diag(nrow(as.matrix(Y)))
            Z <- list(Z=zz, K=Z[[1]])
          }
          if(provided == "Z"){
            #kk <- diag(dim(Z[[1]])[2])
            kk <- diag(dim(Z[[1]])[2])
            attributes(kk)$diagon <- TRUE
            #rownames(kk) <- colnames(Z[[1]]); colnames(kk) <- rownames(kk)
            Z <- list(Z=Z[[1]],K=kk) 
          }
        }else{ # there's 2 matrices in Z
          dido<-lapply(Z, dim) # dimensions of Z and K
          condi<-(dido$Z[2] == dido$K[1] & dido$Z[2] == dido$K[2]) 
          # condition, column size on Z matches with a square matrix K
          if(!condi){
            cat(paste("ERROR! In the",s,"th random effect you have provided or created an incidence \nmatrix with dimensions:",dido$Z[1],"rows and",dido$Z[2],"columns. Therefore the \nvariance-covariance matrix(K) for this random effect expected was a \nsquare matrix with dimensions",dido$Z[2],"x",dido$Z[2]),", but you provided a",dido$K[1],"x",dido$K[2]," matrix \nas a variance-covariance matrix. Please double check your matrices.")
            stop()
          }else{Z=list(Z=Z)}
        }
      }
    }else{
      if(is.null(Z)){ # the user is not using the random part
        cat("Error. No random effects specified in the model. \nPlease use 'lm' or provide a diagonal matrix in Z\ni.e. Zu = list(A=list(Z=diag(length(y))))\n")
        stop()
      }else{
        #stop;
        cat("\nThe parameter 'Z' needs to be provided in a 2-level list structure. \n\nPlease see help typing ?mmer and look at the 'Arguments' section\n")
        cat("\nIf no random effects provided, the model will be fitted using the 'lm' function\n\n")
      }
    }
  }
  
  if(method == "NR"){
    RES <- MNR(Y=Y,X=X,ZETA=Z,R=R,W=W,init=init,iters=iters,tolpar=tolpar,
               tolparinv = tolparinv,draw=draw,silent=silent, 
               constraint = constraint,EIGEND = EIGEND,
               forced=forced,IMP=IMP,restrained=restrained, REML=REML, 
               init.equal = init.equal)
    class(RES)<-c("MMERM")
    #RES$Y.or <- Y
    #RES$ZETA.or <- Z
  }else if(method == "EMMA"){
    if(length(Z)>1){stop("EMMA method only works for one random effect other than error.\n Please select NR or AI methods.", call. = FALSE)}
    RES <- MEMMA(Y=Y,X=X,ZETA=Z,tolpar=tolpar,
                 tolparinv = tolparinv,check.model=check.model,silent=silent)
    class(RES)<-c("MMERM")
  }else if(method == "AI"){
    stop("AI method has been discontinued because of its instability. Try 'NR'.\nSee details in the sommer help page. ", call. = FALSE)
  }else{
    stop("Method not available. See details in the sommer help page.", call. = FALSE)
  }
  
  #   else if(method == "AI"){
  #     RES <- MAI(Y=Y,X=X,ZETA=Z,R=R,init=init,iters=iters,tolpar=tolpar,
  #                tolparinv = tolparinv,draw=draw,silent=silent, 
  #                constraint = constraint,EIGEND = EIGEND,
  #                forced=forced,IMP=IMP,restrained=restrained)
  #     class(RES)<-c("MMERM")
  #   }
  #######$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #######$$$$$$$$$$$$$$$$$$$$$$$$$$$
  layout(matrix(1,1,1))
  return(RES)
}



#### =========== ####
## SUMMARY FUNCTION MMERM #
#### =========== ####
"summary.MMERM" <- function(object, ...) {
  
  #dim(object$u.hat)
  digits = max(3, getOption("digits") - 3)
  #forget <- length(object)
  
  groupss.nn <-do.call("rbind",object$dimos)
  
  colnames(groupss.nn) <- c("Observ","Groups")
  
  if(object$method == "MNR"){
    nano <- object$random.effs
  }else{
    nano <- names(object$ZETA)
  }
  
  if(!is.null(nano)){ # names available
    rownames(groupss.nn) <- nano
    reto <- object$var.comp
    for(i in 1:length(reto)){
      names(reto)[i] <- paste("Var-Covar(",nano[i],")",sep="")
    }
    names(reto)[length(reto)] <- "Var-Covar(Residual)"
  }else{ # no names available
    nano <- paste("u",1:length(object$ZETA),sep=".")
    rownames(groupss.nn) <- nano
    reto <- object$var.comp
    for(i in 1:length(reto)){
      names(reto)[i] <- paste("Var-Covar(",nano[i],")",sep="")
    }
    names(reto)[length(reto)] <- "Var-Covar(Residual)"
  }
  
  LLAIC <- data.frame(as.numeric(object$LL), as.numeric(object$AIC),
                      as.numeric(object$BIC), object$method, object$convergence)
  colnames(LLAIC) = c("logLik","AIC","BIC","Method","Converge")
  rownames(LLAIC) <- "Value"
  
  method=object$method
  #extract fixed effects
  coef <- data.frame((object$beta.hat))#, Std.Error=(matrix(sqrt(diag(object$Var.beta.hat)),ncol=1)), t.value=(matrix((object$beta.hat-0)/sqrt(diag(object$Var.beta.hat)), ncol=1)))
  if(dim(coef)[1] == 1){rownames(coef) <- "Intercept"}
  #if(is.null(rownames(coef))){
  if(!is.null(colnames(object$X))){
    rownames(coef) <- unique(colnames(object$X))
  }
  object$beta.hat
  ## se and t values for fixed effects
  ts <- ncol(coef)
  s2.beta <- diag(as.matrix(object$Var.beta.hat))
  nse.beta <- length(s2.beta)/ts
  inits <- seq(1,length(s2.beta),nse.beta)
  ends <- inits+nse.beta-1
  seti <- list() # stardard errors partitioned by trait
  for(u in 1:ts){
    prox <- data.frame(coef[,u],sqrt(abs(s2.beta[inits[u]:ends[u]])))
    prox$`t value` <- prox[,1]/prox[,2]
    colnames(prox) <- c("Estimate","Std. Error","t value")
    rownames(prox) <- rownames(coef)
    seti[[u]] <- prox
  }
  names(seti) <- colnames(coef)
  
  
  w <- reto#object$var.comp
  
  ########## get the names of the variable combos
  wx <- w[[1]]
  
  if(!is.null(colnames(wx))){
    www <- matrix(NA,dim(wx),dim(wx)) 
    for(u in 1:dim(www)[1]){
      for(uu in 1:dim(www)[1]){
        www[u,uu] <- paste(rownames(wx)[u],rownames(wx)[uu],sep="-")
      }
    }
    bb <- upper.tri(www); diag(bb) <- TRUE
    babasss <- which(bb,arr.ind = TRUE); babasss <- babasss[ order(babasss[,1], babasss[,2]), ]
    sisisi <- www[babasss]
  }
  ########## get the names of the variable combos
  
  if(object$method == "MNR" | object$method == "MAI"){
    #### 111
    ###2222
    zzz <- object$sigma.scaled/sqrt(abs(diag(object$fish.inv)))
    vaross <- lapply(object$var.comp, function(x){
      if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
      babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
      return(x[babas])}else{return(as.matrix(x))}})
    
    if(!is.null(colnames(wx))){
      vaross <- lapply(vaross,function(x){dd <- as.vector(x);names(dd) <- sisisi[1:length(dd)];return(dd)})
      #}
      
      sigma.real <- as.matrix(unlist(vaross))
      
      
      sese <- sigma.real/zzz
      w2 <- data.frame(VarComp=sigma.real, VarCompSE=sese,
                       Zratio=zzz)
    }
  }else{
    ##### 11111
    vaross <- lapply(object$var.comp, function(x){
      if(dim(as.matrix(x))[1]>1){aa <- upper.tri(x); diag(aa) <- TRUE
      babas <- which(aa,arr.ind = TRUE); babas <- babas[ order(babas[,1], babas[,2]), ]
      return(x[babas])}else{return(as.matrix(x))}})
    if(!is.null(colnames(wx))){
      vaross <- lapply(vaross,function(x){dd <- as.vector(x);names(dd) <- sisisi[1:length(dd)];return(dd)})
      #}
      
      sigma.real <- as.matrix(unlist(vaross))
      #sese <- sigma.real/zzz
      w2 <- data.frame(VarComp=sigma.real)#, VarCompSE=sese,Zratio=zzz)
    }
  }
  
  if(!is.null(object$restrained)){
    w2 <- w2[-object$restrained,]
  }
  output <- list(groups=groupss.nn, var.comp=w,var.comp.table=w2, betas=seti, method=method,logo=LLAIC)
  attr(output, "class")<-c("summary.MMERM", "list")
  return(output)
}

"print.summary.MMERM"<-function (x, digits = max(3, getOption("digits") - 3),  ...){
  
  nmaxchar0 <- max(as.vector(unlist(apply(data.frame(rownames(x$var.comp.table)),1,nchar))),na.rm = TRUE)
  
  if(nmaxchar0 < 24){
    nmaxchar0 <- 24
  } # + 26 spaces we have nmaxchar0+26  spaces to put the title
  
  nmaxchar <- nmaxchar0+26 ## add spaces from the 3 columns
  nmaxchar2 <- nmaxchar0+12
  nmaxchar3 <- nmaxchar0+26-46 #round(nmaxchar0/2)
  rlh <- paste(rep("*",round(nmaxchar2/2)),collapse = "")
  rlt <- paste(rep(" ",ceiling(nmaxchar3/2)),collapse = "")
  digits = max(3, getOption("digits") - 3)
  #cat("Information contained in this structure: \n* Results for a multi response model\nDisplayed: \n* Variance-covariance component summaries\nUse the '$' sign to access parameters\n")
  #cat("============================================================")
  cat(paste(rep("=",nmaxchar), collapse = ""))
  #cat("\n   Multivariate Linear Mixed Model fit by REML      \n")
  cat(paste("\n",rlt,"Multivariate Linear Mixed Model fit by REML",rlt,"\n", collapse = ""))
  #cat("***********************  sommer 3.4  ***********************\n")
  cat(paste(rlh," sommer 3.4 ",rlh, "\n", collapse = ""))
  #cat("============================================================")
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\n")
  #print((x$method))
  cat("")
  print(x$logo)#, digits = digits)
  #cat("============================================================")
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nVariance-Covariance components:\n")
  # xx <- x$var.comp.table
  # pxx <- apply(data.frame(rownames(xx)),1,function(x){substr(x,0,33)})
  # dupos <- which(duplicated(pxx))
  # if(length(dupos)>0){
  #   pxx[dupos] <- paste(pxx[dupos],1:length(pxx[dupos]),sep=".")
  # }
  # rownames(xx) <- pxx
  print(x$var.comp.table, digits = digits)
  #cat("============================================================")
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nFixed effects:\n\n")
  print(x$betas, digits = digits)
  #cat("============================================================")
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nGroups and observations:\n")
  print(x$groups, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  #cat("============================================================")
  cat("\nUse the '$' sign to access results and parameters")#\nArguments set to FALSE for multiresponse models:\n'draw', and 'gwas.plots'\n")
}

#### =========== ######
## RESIDUALS FUNCTION #
#### =========== ######
"residuals.MMERM" <- function(object, type="conditional", ...) {
  digits = max(3, getOption("digits") - 3)
  
  if(type=="conditional"){
    output <- object$cond.residuals
    #colnames(output) <- names(object)
  }else{
    output<- output <- object$residuals
    #colnames(output) <- names(object)
  }
  return(output)
}

"print.residuals.MMERM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}
#### =========== ######
## RANEF FUNCTION #
#### =========== ######

"randef" <- function(object) {
  if(class(object)=="MMERM"){
    digits = max(3, getOption("digits") - 3)
    cat("Returning object of class 'list' where each element correspond to one random effect.")
    output <- object$u.hat
  }else{
    stop("Class not recognized.\n",call. = FALSE)
  }
  return(output)
}

#"print.ranef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
#  print((x))
#}

#### =========== ######
## FIXEF FUNCTION #
#### =========== ######


#"fixef.MMERM" <- function(object, ...) {
#  digits = max(3, getOption("digits") - 3)
#  output <- object$beta.hat
#  return(output)
#}

#"print.fixef.MMERM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
#  print((x))
#}

#### =========== ####
## FITTED FUNCTION ##
#### =========== ####

"fitted.MMERM" <- function(object, type="complete", ...) {
  #type="complete" 
  digits = max(3, getOption("digits") - 3)
  if(type=="complete"){
    output<- object$fitted.y
    #colnames(output) <- names(object)
  }else{
    output<- object$fitted.u
    #colnames(output) <- names(object)
  }
  return(output)
}

"print.fitted.MMERM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 


#### =========== ####
## COEF FUNCTION ####
#### =========== ####

"coef.MMERM" <- function(object, ...){
  object$beta.hat
}

"print.coef.MMERM"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 

#### =========== ####
## ANOVA FUNCTION ###
#### =========== ####
anova.MMERM <- function(object, object2=NULL, ...) {
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
      yu <- summary(object)
      yu2 <- summary(object2)
      dis=c(dim(yu$var.comp.table)[1]+dim(object$beta.hat)[1]*dim(object$beta.hat)[2],
            dim(yu2$var.comp.table)[1]+dim(object2$beta.hat)[1]*dim(object2$beta.hat)[2]) # dimensions
      mods=c("mod1","mod2")
      lls=c(object$LL, object2$LL) # likelihoods
      aics=c(object$AIC, object2$AIC) # AIC's
      bics=c(object$BIC, object2$BIC) # AIC's
      vv=which(dis == max(dis))[1] # which has more variance components BIGGER
      vv2=c(1:2)[which(c(1:2)!= vv)] # SMALLER
      LR = (lls[vv] - lls[vv2])
      r.stat= abs(-2*((LR))) # -2(LL1 - LL2)
      df=dis[vv]-dis[vv2]
      chichi=pchisq((r.stat), df, lower.tail=FALSE)
      if(chichi > 1e-5){
        chichi <- round(chichi,5)
      }
      chichi2=paste(as.character(chichi),signifo(chichi), sep=" ")
      ### construct the table
      cat("Likelihood ratio test for mixed models\n")
      cat("==============================================================\n")
      result=data.frame(Df=c(dis[vv],dis[vv2]), AIC=c(aics[vv],aics[vv2]), 
                        BIC=c(bics[vv],bics[vv2]), loLik=c(lls[vv],lls[vv2]), 
                        Chisq=c("",as.character(round(r.stat,5))), 
                        ChiDf=c("",as.character(df)), PrChisq=c("",chichi2 ))
      rownames(result) <- c(mods[vv],mods[vv2])
      print(result)
      cat("==============================================================\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    }
    #}
  }
  #return(result)
}
#### =========== ####
## PLOTING FUNCTION #
#### =========== ####
plot.MMERM <- function(x, stnd=TRUE, ...) {
  digits = max(3, getOption("digits") - 3)
  transp <- function (col, alpha = 0.5){
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,c[3]/255, alpha))
    return(res)
  }
  # std vs residuals, QQplot (std vs teor quantiles), sqrt(std residuals) vs fitted, std res vs leverage = cook's distance
  traits <- ncol(x$fitted.y)
  layout(matrix(1:4,2,2))
  for(i in 1:traits){
    plot(x$fitted.y[x$used.observations,i], scale(x$cond.residuals[,i]), pch=20, col=transp("cadetblue"), ylab="Std Residuals", xlab="Fitted values", main="Residual vs Fitted", bty="n", ...); grid()
    plot(x$fitted.y[x$used.observations,i], sqrt(abs((scale(x$cond.residuals[,i])))), pch=20, col=transp("thistle4"), ylab="Sqrt Abs Std Residuals", xlab="Fitted values", main="Scale-Location",bty="n", ...);grid()
    qqnorm(scale(x$cond.residuals), pch=20, col=transp("tomato1"), ylab="Std Residuals", bty="n",...); grid()
    hat <- x$X%*%solve(t(x$X)%*%x$V.inv%*%x$X)%*%t(x$X)%*%x$V.inv # leverage including variance from random effects H= X(X'V-X)X'V-
    plot(diag(hat), scale(x$cond.residuals), pch=20, col=transp("springgreen3"), ylab="Std Residuals", xlab="Leverage", main="Residual vs Leverage", bty="n", ...); grid()
  }
  #####################
  layout(matrix(1,1,1))
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
  sommer.version = "3.4 (2018-04-01)" # usually 2 months before it expires
  
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
    packageStartupMessage(paste("[]================================================================[]"),appendLF=TRUE)
    packageStartupMessage(paste("[]  Solving Mixed Model Equations in R (sommer) ", sommer.version, "  []",sep=""),appendLF=TRUE)
    packageStartupMessage(paste("[]  ------------ Multivariate Linear Mixed Models --------------  []"),appendLF=TRUE)
    #    packageStartupMessage("[]  ----- Enabling covariance structures in random effects -----  []",appendLF=TRUE)
    packageStartupMessage("[]  Author: Giovanny Covarrubias-Pazaran                          []",appendLF=TRUE)
    packageStartupMessage("[]  Published: PLoS ONE 2016, 11(6):1-15                          []",appendLF=TRUE)
    #    packageStartupMessage("[]  Supported by the Council of Science and Technology (CONACYT)  []", appendLF=TRUE)
    packageStartupMessage("[]  Type 'vignette('sommer.start')' for a short tutorial          []",appendLF=TRUE)
    packageStartupMessage("[]  Type 'citation('sommer')' to know how to cite sommer          []",appendLF=TRUE)
    packageStartupMessage(paste("[]================================================================[]"),appendLF=TRUE)
    packageStartupMessage("UPDATE 'sommer' EVERY 3-MONTHS USING 'install.packages('sommer')'",appendLF=TRUE)
    
    #if(yyy > current){ # yyy < current in CRAN
    #  packageStartupMessage(paste("Version",current,"is now available."),appendLF=TRUE) # version current
    #  packageStartupMessage(paste("Please update 'sommer' installing the new version."),appendLF=TRUE) # version current
    #}
    #print(image(diag(10),main="sommer 3.4"))
  }
  invisible()
}

#.onLoad <- function(library, pkg) {
#  data(x)
#  library(audio)
#  packageStartupMessage(play(x))
#}