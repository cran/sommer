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
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix 
A <- A.mat(CPgeno) # additive relationship matrix
#### look at the data and fit the model
head(CPpheno)
# mix1 <- mmer2(Yield~1,
#               random=~g(id)
#                       + Rowf + Colf, 
#               rcov=~units,
#               G=list(id=A), 
#               data=CPpheno)
# summary(mix1)
#### calculate heritability
# pin(mix1, h1 ~ V1/(V1+V4) )

####=========================================####
#### a more complex example
####=========================================####

# CPpheno$idd <- CPpheno$id; CPpheno$ide <- CPpheno$id;
# head(CPpheno)
# CPgeno[1:4,1:4]

# A <- A.mat(CPgeno) # additive relationship matrix
# D <- D.mat(CPgeno) # dominance relationship matrix
# E <- E.mat(CPgeno) # epistatic relationship matrix
# 
# mix1 <- mmer2(Yield~1,
#               random=~g(id) + g(idd) + g(ide) 
#                       + Rowf + Colf, 
#               rcov=~units,
#               G=list(id=A, idd=D, ide=E), 
#               data=CPpheno)
# summary(mix1)
#### calculate heritability
# pin(mix1, h1 ~ V1/(V1+V6) )
# assumes diag(trait) structure for univariate models

####=========================================####
### if you want to do the same building the matrices
### by your own
### (not recommended unless mmer2 can't do your model)
####=========================================####

# data(CPdata)
# #### look at the data
# head(CPpheno)
# CPgeno[1:5,1:5]
# #### fit a model including additive and dominance effects
# y <- CPpheno$color
# Za <- diag(length(y))
# Zd <- diag(length(y))
# Ze <- diag(length(y))
# A <- A.mat(CPgeno) # additive relationship matrix
# D <- D.mat(CPgeno) # dominance relationship matrix
# E <- E.mat(CPgeno) # epistatic relationship matrix
# 
# y.trn <- y # for prediction accuracy
# ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
# y.trn[ww] <- NA
# ####================####
# #### ADDITIVE MODEL ####
# ####================####
# ETA.A <- list(add=list(Z=Za,K=A))
# ans.A <- mmer(Y=y.trn, Z=ETA.A)
# cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
# ####=========================####
# #### ADDITIVE-DOMINANCE MODEL ####
# ####=========================####
# ETA.AD <- list(add=list(Z=Za,K=A), dom=list(Z=Zd,K=D))
# ans.AD <- mmer(Y=y.trn, Z=ETA.AD)
# cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
# ### greater accuracy !!!! 4 percent increment!!
# ### we run 100 iterations, 4 percent increment in general
# ####===================================####
# #### ADDITIVE-DOMINANCE-EPISTASIS MODEL ####
# ####===================================####
# ETA.ADE <- list(add=list(Z=Za,K=A), dom=list(Z=Zd,K=D), epi=list(Z=Ze,K=E))
# ans.ADE <- mmer(Y=y.trn, Z=ETA.ADE)
# cor(ans.ADE$fitted.y[ww], y[ww], use="pairwise.complete.obs")
# 
# summary(ans.A)
# summary(ans.AD)
# summary(ans.ADE)
# #### adding more effects doesn't necessarily increase prediction accuracy!

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
#### Univariate vs Multivariate models
####=========================================####
####=========================================####

# data(CPdata)
# head(CPpheno)
# CPgeno[1:4,1:4]
# #### create the variance-covariance matrix 
# A <- A.mat(CPgeno)
# 
# CPpheno2 <- CPpheno
# ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
# CPpheno2[ww,"color"] <- NA
# 
# ####==================####
# #### univariate model ####
# ####==================####
# ans.u <- mmer2(color~1,
#               random=~ g(id) 
#               + Rowf + Colf, 
#               rcov=~units,
#               G=list(id=A), 
#               data=CPpheno2)
# cor(ans.u$u.hat$`g(id)`[ww,], CPpheno[ww,"color"], use="pairwise.complete.obs")
# ####====================####
# #### multivariate model ####
# ####     2 traits       ####
# ####====================####
# #### be patient take some time
# ans.m <- mmer2(cbind(Yield,color)~1,
#                random=~ us(trait):g(id) 
#                + diag(trait):Rowf + diag(trait):Colf, 
#                rcov=~ us(trait):units,
#                G=list(id=A), 
#                data=CPpheno2)
# cor(ans.m$u.hat$`g(id)`[ww,2], CPpheno[ww,"color"], use="pairwise.complete.obs")


####=========================================####
### if you want to do the same building the matrices
### by your own
### (not recommended unless mmer2 can't do your model)
####=========================================####

# data(CPdata)
# ### look at the data
# head(CPpheno)
# CPgeno[1:5,1:5]
# ## fit a model including additive and dominance effects
# Za <- diag(dim(CPpheno)[1])
# A <- A.mat(CPgeno) # additive relationship matrix
# 
# CPpheno2 <- CPpheno
# ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
# CPpheno2[ww,"color"] <- NA
# ####==================####
# #### univariate model ####
# ####==================####
# ETA.A <- list(add=list(Z=Za,K=A))
# ans.u <- mmer(Y=CPpheno2$color, Z=ETA.A)
# cor(ans.u$u.hat$`g(id)`[ww,], CPpheno[ww,"color"], use="pairwise.complete.obs")
# ####====================####
# #### multivariate model ####
# ####     2 traits       #### 
# ####====================####
# #### be patient take some time
# ETA.B <- list(add=list(Z=Za,K=A))
# ans.m <- mmer(Y=CPpheno2[,c("Yield","color")], Z=ETA.B)
# cor(ans.m$u.hat$`g(id)`[ww,2], CPpheno[ww,"color"], use="pairwise.complete.obs")

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }
