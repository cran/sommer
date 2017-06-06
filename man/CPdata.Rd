\name{CPdata}
\alias{CPdata}
\docType{data}
\title{
Genotypic and Phenotypic data for a CP population
}
\description{
A CP population or F1 cross is the designation for a cross between 2 highly heterozygote individuals; i.e. humans, fruit crops, bredding populations in recurrent selection.

This dataset contains phenotpic data for 363 siblings for an F1 cross. These are averages over 2 environments evaluated for 4 traits; color, yield, fruit average weight, and firmness. The columns in the CPgeno file are the markers whereas the rows are the individuals. The CPpheno data frame contains the measurements for the 363 siblings, and as mentioned before are averages over 2 environments.

}
\usage{data("CPdata")}
\format{
  The format is:
 chr "CPdata"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
 
}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(CPdata)
CPpheno <- CPdata$pheno[,-c(1:4)]
CPgeno <- CPdata$geno
#### look at the data
head(CPpheno)
CPgeno[1:5,1:5]
#### fit a model including additive and dominance effects
y <- CPpheno$color
Za <- diag(length(y))
Zd <- diag(length(y))
Ze <- diag(length(y))
A <- A.mat(CPgeno) # additive relationship matrix
#D <- D.mat(CPgeno) # dominant relationship matrix
#E <- E.mat(CPgeno) # epistatic relationship matrix

y.trn <- y # for prediction accuracy
ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
y.trn[ww] <- NA
####================####
#### ADDITIVE MODEL ####
####================####
#ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- mmer(Y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================####
#### ADDITIVE-DOMINANCE MODEL ####
####=========================####
#ETA.AD <- list(add=list(Z=Za,K=A), dom=list(Z=Zd,K=D))
#ans.AD <- mmer(Y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
### greater accuracy !!!! 4 percent increment!!
### we run 100 iterations, 4 percent increment in general
####===================================####
#### ADDITIVE-DOMINANCE-EPISTASIS MODEL ####
####===================================####
#ETA.ADE <- list(add=list(Z=Za,K=A), dom=list(Z=Zd,K=D), epi=list(Z=Ze,K=E))
#ans.ADE <- mmer(Y=y.trn, Z=ETA.ADE)
#cor(ans.ADE$fitted.y[ww], y[ww], use="pairwise.complete.obs")

#summary(ans.A)
#summary(ans.AD)
#summary(ans.ADE)
#### adding more effects doesn't necessarily increase prediction accuracy!

####=========================================####
#### RUN AS GWAS MODEL         
####=========================================####
#ans.GWAS <- mmer(Y=y, Z=ETA.A, W=CPgeno)
#### or if you have a map
#my.map <- CPdata$map
#ans.GWAS <- mmer(Y=y, Z=ETA.A, W=CPgeno, map=my.map)


############################################
############################################
############################################
############################################
############################################
############################################


####=========================================####
####=========================================####
#### EXAMPLE 2
#### Genomic prediction
#### using multiple responses in parallel
####=========================================####
####=========================================####
#data(CPdata)
#CPpheno <- CPdata$pheno[,-c(1:4)]
#CPgeno <- CPdata$geno
### look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]

#### fit a model including additive effects
#y <- CPpheno$color
#Za <- diag(length(y))
#A <- A.mat(CPgeno) # additive relationship matrix

#y.trn <- CPpheno # for prediction accuracy
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#y.trn[ww,] <- NA
#ETA.A <- list(add=list(Z=Za,K=A))
#********************************************************
#### WINDOWS
#ans.A <- mmer(Y=y.trn, Z=ETA.A, method="EMMA")
#### MAC
#ans.A <- mmer(Y=y.trn, Z=ETA.A, method="EMMA",n.cores=4)
#********************************************************

### 1) Extract the fitted values for the simulated data
#fitt <- as.matrix(do.call("cbind",lapply(ans.A, function(x,w){x$fitted.y[w]},w=ww)))
### 2) See the correlation with the original data masked
#res <- apply(data.frame(1:4),1,function(x,y,z){
#  cor(y[,x],z[,x],use="pairwise.complete.obs")},
#             y=fitt,z=CPpheno[ww,])
#res

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 3
#### Genomic prediction
#### Univariate vs Multivariate models
####=========================================####
####=========================================####
#data(CPdata)
#CPpheno <- CPdata$pheno[,-c(1:4)]
#CPgeno <- CPdata$geno
### look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
## fit a model including additive and dominance effects
#Za <- diag(dim(CPpheno)[1])
#A <- A.mat(CPgeno) # additive relationship matrix

#CPpheno2 <- CPpheno
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#CPpheno2[ww,"color"] <- NA
####==================####
#### univariate model ####
####==================####
#ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- mmer(Y=CPpheno2$color, Z=ETA.A)
#cor(ans.A$fitted.y[ww], CPpheno[ww,"color"], use="pairwise.complete.obs")
####====================####
#### multivariate model ####
####     4 traits       #### 
####====================####
#### be patient take some time
#ETA.B <- list(add=list(Z=Za,K=A))
#ans.B <- mmer(Y=CPpheno2, Z=ETA.B, MVM=TRUE)
#cor(ans.B$fitted.y[ww,"color"], CPpheno[ww,"color"], use="pairwise.complete.obs")


############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 4
#### Cross validation example
####=========================================####
####=========================================####
## get data
# data(CPdata)
# CPpheno <- CPdata$pheno[,-c(1:4)]
# CPgeno <- CPdata$geno
# head(CPpheno)
# CPgeno[1:5,1:5]
# 
# y <- CPpheno$color
# Za <- diag(length(y))
# A <- A.mat(CPgeno)
# 
# ## for 5-fold cross validation
# rownames(CPpheno) <- rownames(CPgeno)
# 
# cv.cors <- list()
# for(j in 1:50){ # 50 rounds of 5-fold CV
#   
#   base <- rownames(CPgeno) # base names of the genotypes
#   ## define the cross validation groups
#   CV.scheme <- list()
#   for(i in 1:5){
#     if(i ==5){totake <- length(base)}else{totake<-round(dim(Za)[1]/5)}
#     CV.scheme[[i]] <- sample(base,totake) # delete data for 1/5 of the population
#     base <- setdiff(base,CV.scheme[[i]])
#   }
#   ## apply the model to all groups
#   coros <- lapply(CV.scheme,function(x){
#     y.trn<- CPpheno
#     y.trn[x,"Yield"] <- NA
#     ETA.A <- list(add=list(Z=Za,K=A))
#     ans.A <- mmer(Y=y.trn[,"Yield"], Z=ETA.A, method="EMMA", silent = TRUE)
#     blup <- ans.A$u.hat$add;rownames(blup) <- rownames(CPpheno)
#     cor(blup[x,1], CPpheno[x,"Yield"], use="pairwise.complete.obs")
#   })
#   print(j)
#   cv.cors[[j]] <- unlist(coros)
# }
# 
# boxplot(unlist(cv.cors),ylim=c(0,1))
# 
# ans.A <- mmer(Y=CPpheno[,"Yield"], Z=ETA.A, method="EMMA", silent = TRUE)
# (h2<-ans.A$var.comp[1,1]/sum(ans.A$var.comp))
# sqrt(h2)

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }
