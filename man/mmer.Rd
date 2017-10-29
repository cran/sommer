\name{mmer}
\alias{mmer}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations in \strong{R}

}
\description{
Fits a multivariate or univariate linear mixed model by likelihood methods (REML). It has been \strong{implemented to work directly with incidence and variance covariance matrices} for each random effect. For a formula-based specification please use the \code{\link{mmer2}} function. The optimization methods are; Direct-Inversion Newton-Raphson (NR) (Tunnicliffe 1989; Lee et al. 2015), and Efficient Mixed Model Association (EMMA) (Kang et al. 2008). These algorithms are \strong{intended to be used for problems of the type p > n}. For problems with sparse covariance structures, or problems of the type n > p, the MME-based algorithms are faster and we recommend to shift to the use of MME-based software (i.e. lme4, breedR, or asreml-R). 

The package also provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices to model the covariance among genotypes in plant and animal breeding problems.

\strong{For a short tutorial} of how to perform different genetic analysis with sommer please type:

\strong{vignette("sommer.start")} or \strong{vignette("sommer")}

}
\usage{
mmer(Y,X=NULL,Z=NULL,R=NULL,method="NR",init=NULL,iters=20,
    tolpar=1e-3, tolparinv=1e-6,draw=FALSE,silent=FALSE, 
    constraint=TRUE, EIGEND=FALSE, forced=NULL, IMP=FALSE, 
    complete=TRUE, check.model=TRUE, restrained=NULL, REML=TRUE)
}
\arguments{
  \item{Y}{a matrix or data frame of response variables}
  \item{X}{an incidence matrix for fixed effects.}
  \item{Z}{incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. \strong{This needs to be provided as a 2-level list structure}. For example:
  
  '
  
  ETA <- list(
  
  A=list(Z=Z1, K=K1),
  
  B=list(Z=Z2, K=K2),
  
  C=list(Z=Z3, K=K3)
  
  ) 
  
  
  '
  
  makes a 2 level list for 3 the random effects A, B and C, stored in a variable we call ETA. The general idea is that each random effect is a list, i.e. A=list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix for the random effect, \strong{if K is not provided is assumed an identity matrix} conferring independence. 
  
  '
  
  PLEASE remember to \strong{use the names Z and K for all random effects} when you provide your matrices, \strong{that's the only way the program distinguishes between a Z or a K matrix}. 
  
  '
  
  To provide extra detail, I'll rephrase it; when moving to situations of more than one random effect, you need to build a list for each random effect, and at the end everything gets joined in a list as well (BGLR type of format). Is called a 2-level list, i.e. A=list(Z=Z1, K=K1) and B=list(Z=Z2, K=K2) refers to 2 random effects and they should be put together in a list:
  
  '
  
  ETA <- list( A=list(Z=Z1, K=K1), B=list(Z=Z1, K=K1) )
  
  '
  
  Now you can fit your model as:
  
  '
  
  mod1 <- mmer(Y=y, Z=ETA) 
  
  '
  
  You can see the examples at the bottom to have a clearer idea how to fit your models. }
 \item{R}{list of R matrices for residual effects. If multiple residual structures are provided a direct sum will be applied and the user has to keep in mind that all matrices have to be of the order NxN (neing N the number of observations) even though they will be full of zeros in certain parts of the matrix. Good knowledge of linear algebra is required to use this argument.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. Sommer now focus only in 2 algorithms; Direct-inversion Newton-Raphson \strong{NR} (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2015), and \strong{EMMA} efficient mixed model association (Kang et al. 2008). Average Information has been discontinued because of the instability for very complex models requiring EM-type of algorithms to work.}
  \item{init}{initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional. It has to be provided as a list, where each list element is one variance component and if multitrait model is pursued each element of the list is a matrix of variance covariance components among traits.}
  \item{iters}{Maximum number of iterations allowed.  Default value is
    15.}
  \item{tolpar}{Convergence criteria.}
  \item{tolparinv}{tolerance parameter for matrix inverse}
    \item{draw}{a TRUE/FALSE value indicating if a plot showing the progress of the variance components should be drawn or not. Default is FALSE}
    \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar as iterations progress. The default is FALSE, which means is not silent and will display the progress bar.}
    \item{constraint}{a TRUE/FALSE value indicating if the function should apply the boundary constraint indicating that variance components that are zero should be restrained.}
    \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this method are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM. If you use this feature keep in mind that residuals and fitted values cannot be used given the roation of the data.}
    \item{forced}{a list of values for variance-covariance components to be used if the user wants to force those values. The format is the same than the init argument.}
    \item{IMP}{a TRUE/FALSE statement if the function should impute the Y matrix/dataframe with the median value or should get rid of missing values for the estimation of variance components. Default is FALSE (although please see the 'complete' argument.)}
   \item{complete}{a TRUE/FALSE statement to indicate if the function should impute the cases where at least for one trait there's an observation. Default is TRUE in order to use all data points where at least one observation exist across traits.}
   \item{check.model}{a TRUE/FALSE statement to indicate if the function should check the input parameters from the user.}
   \item{restrained}{a numeric argument specifying which variance-covariance parameters should be restrained. Please see \code{\link{vctable}} and \code{\link{vctable.help}} functions if you want to use it although you shouldn't since is made for the \code{\link{mmer2}} function.}
   \item{REML}{a TRUE/FALSE value to indicate if REML or ML should be used for optimization. Not functional yet. Only REML available.}
    
}
\details{

The package has been equiped with several datasets to learn how to use the sommer package: 

* \code{\link{HDdata}} and \code{\link{FDdata}} datasets have examples to fit half and full diallel designs. 

* \code{\link{h2}} to calculate heritability

* \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets to perform genomic prediction in hybrid single crosses 

* \code{\link{wheatLines}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link{CPdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. 
* \code{\link{PolyData}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link{gryphondata}} data contains an example of an animal model including pedigree information.

* \code{\link{BTdata}} dataset contains an animal (birds) model.

If you need to use pedigree you need to convert your pedigree into a relationship matrix (use the `getA` function from the pedigreemm package). Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the mmer-type of functions. Additional functions for genetic analysis have been included such as False Discovery Rate calculation (\code{\link{fdr}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}), phasing F1 or CP populations (\code{\link{phase.F1}}). We have also included the famous \code{\link{transp}} function to get transparent colors. Useful functions for analyzing such designs are included such as the \code{\link{blocker}} and \code{\link{fill.design}}.

\strong{MODELS ENABLED}

For details about the models enabled please check the \code{\link{sommer}} help page.

}
\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{var.comp}{a list with the values of the variance-covariance components with one list element for each random effect.}
\item{V.inv}{the inverse of the phenotypic variance matrix V^- = (ZGZ+R)^-1}
\item{u.hat}{a list (one element for each random effect) with a data frame for trait BLUPs.}
\item{Var.u.hat}{a list (one element for each random effect) with the variance-covariance matrix for trait BLUPs.}
\item{PEV.u.hat}{a list (one element for each random effect) with the predicted error variance matrix for trait BLUPs.}
\item{beta.hat}{a data frame for trait BLUEs (fixed effects).}
\item{Var.beta.hat}{a variance-covariance matrix for trait BLUEs}
\item{fish.inv}{inverse of the Fisher's information.}
\item{residuals}{Residual values e = Y - XB}
\item{cond.residuals}{Conditional residual e.cond = Y - (XB + ZU)}
\item{LL}{LogLikelihood}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{X}{ incidence matrix for fixed effects}
\item{dimos}{dimnensions for incidence matrix for random effects}
\item{sigma.scaled}{scaled variance covariance components}
\item{fitted.y}{ Fitted values y.hat=XB+Zu}
\item{fitted.u}{ Fitted values only across random effects u.hat=Zu.1+....+Zu.i}
\item{ZETA}{Original incidence and variance covariance matrices used in the model fit.}
\item{K}{ variance-covariance matrix for random effects. If more than one random effect this will be the diagonal binding of individual K matrices.}
\item{method}{The method for extimation of variance components specified by the user.}
\item{forced}{a TRUE/FALSE statement indicating if user used the forced argument.}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{restrained}{table of restrained parameters.}

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744 

Covarrubias-Pazaran G et al. Validating multivariate genomic selection and genome wide association in American cranberry. Submitted to the Plant Genome (2016).

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42:355-360.

}
\author{
Giovanny Covarrubias-Pazaran
}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####

####=========================================####
####=========================================####
#### EXAMPLE 1
#### breeding values with 1 variance component
####=========================================####
####=========================================####

####=========================================####
#### simulate genotypic data
#### random population of 200 lines with 1000 markers
####=========================================####
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
####=========================================####
#### simulate phenotypes
####=========================================####
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
M <- t(M)
####=========================================####
#### fit the model
####=========================================####
Z1 <- diag(length(y))
ETA <- list( add=list(Z=Z1, K=A.mat(M)) )
ans <- mmer(Y=y, Z=ETA)
summary(ans)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 2
#### breeding values with 3 variance components
#### Hybrid prediction
####=========================================####
####=========================================####
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K
y <- hybrid2$Yield
X1 <- model.matrix(~ Location, data = hybrid2);dim(X1)
Z1 <- model.matrix(~ GCA1 -1, data = hybrid2);dim(Z1)
Z2 <- model.matrix(~ GCA2 -1, data = hybrid2);dim(Z2)
Z3 <- model.matrix(~ SCA -1, data = hybrid2);dim(Z3)

colnames(Z1) <- levels(hybrid2$GCA1)
colnames(Z2) <- levels(hybrid2$GCA2)
colnames(Z3) <- levels(hybrid2$SCA)
####=========================================####
#### Realized IBS relationships for set of parents 1
####=========================================####
#K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1) 
####=========================================####
#### Realized IBS relationships for set of parents 2
####=========================================####
#K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)
####=========================================####
#### Realized IBS relationships for cross 
#### (as the Kronecker product of K1 and K2)
####=========================================####
#S <- kronecker(K1, K2) ; dim(S)   
#rownames(S) <- colnames(S) <- levels(hybrid2$SCA)

#ETA <- list(GCA1=list(Z=Z1, K=K1), 
#            GCA2=list(Z=Z2, K=K2), 
#            SCA=list(Z=Z3, K=S)
#            )
#ans <- mmer(Y=y, X=X1, Z=ETA)
#ans$var.comp
#summary(ans)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 3
#### use data from example 2
#### COMPARE WITH MCMCglmm
####=========================================####
####=========================================####

####=========================================####
#### the same model run in MCMCglmm:
####=========================================####
#library(MCMCglmm)
# pro <- list(GCA1 = as(solve(K1), "sparseMatrix"), GCA2 = as(solve(K2),
#      + "sparseMatrix"), SCA = as(solve(S), "sparseMatrix") )
#system.time(mox <- MCMCglmm(Yield ~ Location, random = ~ GCA1 + GCA2 + SCA,
#      + data = hybrid2, verbose = T, ginverse=pro))
## Takes 7:13 minutes in MCMCglmm, in sommer only takes 3.4 seconds

####=========================================####
#### it is also possible to do GWAS for hybrids, separatting 
#### and accounting for effects of GCA1, GCA2, SCA
####=========================================####

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 4
#### use data from example 2
#### COMPARE WITH cpgen
####=========================================####
####=========================================####

#Z_list = list(Z1,Z2,Z3)
#G_list = list(solve(K1), solve(K2), solve(S))
#fit <- clmm(y = y, Z = Z_list, ginverse=G_list, niter=15000, burnin=5000)
####=========================================####
#### inspect results and notice that variance 
#### components were NOT estimated correctly!!
#### also takes longer and no user-friendly 
####=========================================####
#str(fit)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 5
#### COMPARE WITH pedigreemm example
####=========================================####
####=========================================####

#library(pedigreemm)
#A <- as.matrix(getA(pedCowsR))
####=========================================####
##### or using mmer2 would look:
####=========================================####
#head(milk)
#fm3 <- mmer2(fixed=milk ~ 1, random = ~ g(id) + herd, 
#             G=list(id=K1), data=milk)
#summary(fm3)
####=========================================####
#### Try pedigreemm but takes longer, 
#### is an extension of lme4
####=========================================####
#fm2 <- pedigreemm(milk ~ (1 | id) + (1 | herd),data = milk, pedigree = list(id= pedCowsR))
#plot(fm3$u.hat$`g(id)`[,1], ranef(fm2)$id[,1])
#plot(fm3$u.hat$herd[,1], ranef(fm2)$herd[,1])
####=========================================####
#### a big data frame with 3397 rows and 1359 animals analyzed
####=========================================####

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 6
#### PREDICTING SPECIFIC PERFORMANCE 
#### within biparental population    
####=========================================####
####=========================================####

#data(CPdata)
## look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
####=========================================####
#### fit a model including additive and dominance effects
####=========================================####
#y <- CPpheno$color
#Za <- diag(length(y))
#Zd <- diag(length(y))
#A <- A.mat(CPgeno)
#D <- D.mat(CPgeno)

#y.trn <- y # for prediction accuracy
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#y.trn[ww] <- NA

####================####
#### ADDITIVE MODEL ####
####================####
#ETA.A <- list(list(Z=Za,K=A))
#ans.A <- mmer(Y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================####
#### ADDITIVE-DOMINANT MODEL ####
####=========================####
#ETA.AD <- list(list(Z=Za,K=A),list(Z=Zd,K=D))
#ans.AD <- mmer(Y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
### greater accuracy !!!! 4 percent increment!!
### we run 100 iterations, 4 percent increment in general
####===================================####
#### ADDITIVE-DOMINANT-EPISTATIC MODEL ####
####===================================####
#ETA.ADE <- list(list(Z=Za,K=A),list(Z=Zd,K=D),list(Z=Ze,K=E))
#ans.ADE <- mmer(Y=y.trn, Z=ETA.ADE)
#cor(ans.ADE$fitted.y[ww], y[ww], use="pairwise.complete.obs")
#### adding more effects doesn't necessarily increase prediction accuracy!

########## NOTE
## nesting in R is indicated as 
## assume blocks nested in locations
## Loc + Block/Loc
## is the same than
## Loc + Block + Loc:Block

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 7
#### MULTIVARIATE MODEL & GENETIC CORRELATION
####=========================================####
####=========================================####

# data(CPdata)
# ### look at the data
# head(CPpheno)
# CPgeno[1:5,1:5]
# ### fit a model including additive and dominance effects
# Y <- CPpheno[,c("color","Yield")]
# Za <- diag(dim(Y)[1])
# A <- A.mat(CPgeno) # additive relationship matrix
# ###================####
# ### ADDITIVE MODEL ####
# ###================####
# ETA.A <- list(add=list(Z=Za,K=A))
# ans.A <- mmer(Y=Y, Z=ETA.A)
# summary(ans.A)
# 
# ### genetic variance covariance
# gvc <- ans.A$var.comp$add
# ### extract variances (diagonals) and get standard deviations
# sd.gvc <- as.matrix(sqrt(diag(gvc))) 
# ### get possible products sd(Vgi) * sd(Vgi')
# prod.sd <- sd.gvc %*% t(sd.gvc) ###remove the \
# ### genetic correlations cov(gi,gi')/[sd(Vgi) * sd(Vgi')]
# (gen.cor <- gvc/prod.sd)
# ### heritabilities
# (h2 <- diag(gvc) / diag(cov(Y, use = "complete.obs")))


############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 8
#### Genomic prediction
#### Univariate vs Multivariate models
####=========================================####
####=========================================####
#data(CPdata)
#### look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
#### fit a model including additive and dominance effects
#Za <- diag(dim(CPpheno)[1])
#A <- A.mat(CPgeno) # additive relationship matrix

#CPpheno2 <- CPpheno
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#CPpheno2[ww,"Yield"] <- NA

# ###==================####
# ### univariate model ####
# ###==================####
# ETA.A <- list(add=list(Z=Za,K=A))
# ans.A <- mmer(Y=CPpheno2$Yield, Z=ETA.A, method="NR")
# ans.A$var.comp
# cor(ans.A$fitted.y[ww], CPpheno[ww,"Yield"], use="pairwise.complete.obs")
# ###====================####
# ### multivariate model ####
# ###====================####
# ans.B <- mmer(Y=CPpheno2[,c("Yield","color","FruitAver")], Z=ETA.A)
# ans.B$var.comp
# cor(ans.B$fitted.y[ww,"Yield"], CPpheno[ww,"Yield"], use="pairwise.complete.obs")

}