\name{mmer}
\alias{mmer}
\title{Mixed Model Equations in R}
\description{
This function is the core of the package and solves the mixed model equations proposed by Henderson (1975). It has been implemented to work with incidence matrices and variance covariance matrices for each random effect. In the details we will explain the methods implemented by this function. Currently 3 methods are supported; "EMMA" efficient mixed model association (Kang et al. 2008), "AI" average information (Gilmour et al. 1995; Lee et al. 2015), and "EM" expectation maximization (Searle 1993; Bernardo 2010). The EMMA method is implemented when only one variance component other than the error variance component (Var(e)) is estimated, is based on optimizing the likelihood function (see details). On the other hand when more than one variance component needs to be estimated the "AI" and "EM"" methods should be used. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominant (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices that have been shown to increase prediction accuracy. The package provides flexibility to fit other genetic models such as full and half diallel models as well, see \code{\link{hdm}} function information to see how to fit those models using sommer. 
}
\usage{
mmer(y, X=NULL, Z=NULL, W=NULL, R=NULL, method="AI", REML=TRUE, 
     iters=50, draw=FALSE, init=NULL, n.PC=0, P3D=TRUE, 
     min.MAF = 0.05, silent=FALSE, family=NULL)
}
\arguments{
  \item{y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects related to environmental effects. This has to be provided as a matrix, NOT in a list structure.}
  \item{Z}{incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. THIS NEEDS TO BE PROVIDED AS A 2-LEVEL LIST STRUCTURE. For example:
  
  .
  
  Z=list(list(Z=Z1, K=K1),list(Z=Z2, K=K2),list(Z=Z3, K=K3)) 
  
  .
  
  makes a 2 level list for 3 random effects. The general idea is that each random effect is a list, i.e. list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix for the random effect, if K is not provided is assumed an incidence matrix. 
  
  .
  
  PLEASE remember to use the names Z and K FOR ANY RANDON EFFECT when you provide your matrices since that's the only way the program distinguishes between a Z or a K matrix. 
  
  .
  
  To provide extra detail, rephrasing; when moving to more than one random effect we need to make a list for each random effect, and at the end everything goes in a list as well (BGLR type of format). We call it a 2-level list, i.e. list(Z=Z1, K=K1) and list(Z=Z2, K=K2) refers to 2 random effects and they should be put together in a list:
  
  ETA = list( list(Z=Z1, K=K1), list(Z=Z1, K=K1) )
  
  Now you can fit your model as:
  
  mod1 <- mmer(y=y, Z=ETA) 
  
  You can see the examples at the bottom to have a clearer idea how to fit your models.
  
  }
  \item{W}{an incidence matrix for extra fixed effects and only to be used if GWAS will be performed and markers effects will be estimated as fixed effects according to Yu et al. (2006). Theorically X and W are both fixed effects, but they are separated to perform GWAS easier and allow the program to recognize the markers from other fixed factors such as environmental factors. This has to be provided as a matrix same than X.}
  \item{R}{a matrix for variance-covariance structures for the residuals, i.e. for longitudinal data. if not passed is assumed an identity matrix.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. The package currently is supported by 3 algorithms; "EMMA" efficient mixed model association (Kang et al. 2008), "AI" average information (Gilmour et al. 1995; Lee et al. 2015), and "EM" expectation maximization (Searle 1993; Bernardo 2010). The default method is average information "AI" because of its ability to handle multiple random effects and its greater speed compared to "EM" which can handle multiple random effects but it is much slower and is not based on REML. On the other hand the "EMMA" method has been implemented but the user needs to keep in mind that such algorithm is based on the theory of ridge regression and the ridge parameter allows the estimation of only a single variance component other than the error.}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM algorithm is selected. There is no rule of tumb for the number of iterations. The default value is 50 iterations or EM steps, but could take less or much longer than that.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is FALSE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
  \item{init}{an vector of initial values for the EM algorithm if this is the method selected. If not provided the program uses a starting values the variance(y)/#random.eff which is usually a good starting value.}
  \item{n.PC}{Number of principal components to include as fixed effects. Default is 0 (equals K model).}
  \item{P3D}{When P3D=TRUE, variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}
  \item{min.MAF}{a scalar value between 0-1 indicating what is theminor allele frequency to be allowed for a marker during a GWAS analysis when providing the matrix W of markers. In general is known that results for markers with alleles with MAF < 0.05 are not reliable unless sample size is big enough.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar and poems (see \code{\link{poe}} function) while working or should not be displayed. The default is FALSE, which means is not silent and will display the progress bar and a short poem to help the scientist (and me haha) remember that life is more than analyzing data.}
  \item{family}{a family object to specify the distribution of the response variable. The program will only use the link function to transform the response. For details see \code{\link{family}} help page. The argument would look something like this; family=poisson(), or family=Gamma(), etc. For more sophisticated models please look at lme4 package from Douglas Bates. NOT IMPLEMENTED YET.}

}
\details{
The package has been developed to provide R users with code to understand how most common algorithms in mixed model analysis work related to genetics field, but also allowing to perform their real analysis. This package allows the user to calculate the variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects. This program focuses in the mixed model of the form:

.

y = Xb + Zu + e ..............where Zu can contain several random effects

.

b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, K*sigma(u)]  .......where: K*sigma(u) = G

e ~ N[0, I*sigma(e)]  .......where: I*sigma(e) = R

y ~ N[Xb, var(Zu+e)]   ......where;

var(y) = var(Zu+e) = ZGZ+R = V which is the phenotypic variance

.

The function allows the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with the following design:

.

fixed = only intercept...........................b ~ N[b.hat, 0]

random = GCA1 + GCA2 + SCA.................u ~ N[0, G]       

.

where G is:

.

|K*sigma2(gca1).....................0..........................0.........|
    
|.............0.............S*sigma2(gca2).....................0.........| = G

|.............0....................0......................W*sigma2(sca)..|

.

This would be specified in the function as:

.

X1 <- matrix(1,length(y),1)      incidence matrix for intercept only

ETA1 <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=K3))      for 3 random effects

.

where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively and K1, K2, K3 their var-cov matrices. Now the fitted model will look like:

.

ans <- mmer(y=y, X=X1, Z=ETA1)

.

FOR DETAILS ON HOW THE "AI", EM" AND "EMMA" ALGORITHMS WORK PLEASE REFER TO \code{\link{AI}} , \code{\link{EM}} AND \code{\link{EMMA}}

In addition, the package contains a very nice function to plot genetic maps with numeric variable or traits next to the LGs, see the \code{\link{map.plot2}} function to see how easy can be done. The package contains other functions:

\code{\link{transp}} function transform a vector of colors in transparent colors.

\code{\link{fdr}} calculates the false discovery rate for a vector of p-values.

\code{\link{A.mat}} is a wrapper of the A.mat function from the rrBLUP package.

\code{\link{D.mat}} calculates the dominant relationship matrix.

\code{\link{E.mat}} calculates de epistatic relationship matrix.

\code{\link{score.calc}} is a function that can be used to calculate a -log10 p-value for a vector of BLUEs for marker effects.

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{ranef}}, \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using this function which is the core of the sommer package.

}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{

\item{$var.com}{ a vector with the values of the variance components estimated}
\item{$V.inv}{ a matrix with the inverse of the phenotypic variance V = ZGZ+R, V^-1}
\item{$u.hat}{ a vector with BLUPs for random effects}
\item{$Var.u.hat}{ a vector with variances for BLUPs}
\item{$PEV.u.hat}{ a vector with predicted error variance for BLUPs}
\item{$beta.hat}{ a vector for BLUEs of fixed effects}
\item{$Var.beta.hat}{ a vector with variances for BLUEs}
\item{$X}{ incidence matrix for fixed effects}
\item{$score}{ the -log10(p-value) for each marker if a GWAS model is fitted by specifying the W parameter in the model.}

In addition, we have included a couple of random poems from Latin American writers to help the scientist (an me haha) remember from time to time that life is more than analyzing data. You can always silence this feature by setting the argument silent=TRUE, which will avoid the program to display the poems. If you want to contribute with a poem, phrase or short citation for future versions of sommer, feel free to send it to me to: 

covarrubiasp@wisc.edu

Please share your ideas and code, future generations of scientists can be better if we are not greedy sharing our knowledge. Feel free to use my code for your own software! good luck with your analysis.

}
}
\references{

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

}
\author{
Giovanny Covarrubias-Pazaran
}
\examples{
############################################
############################################
# breeding values with 1 variance component
############################################
############################################

##== simulate data
##== random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
##== phenotypes
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
M <- t(M)
##== fit the model
Z1 <- diag(length(y))
ETA <- list( list(Z=Z1, K=A.mat(M)))
ans <- mmer(y=y, Z=ETA, method="EMMA")

############################################
############################################
# GWAS with 1 variance component and one A matrix
############################################
############################################

ETA <- list( list(Z=Z1, K=A.mat(M))) # random effects for genotypes
ETA2 <- M # markers as fixed effects 
# RUN IT:
#ans <- mmer(y=y, Z=ETA, W=ETA2, method="EMMA")

############################################
############################################
# breeding values with 3 variance components
# hybrid prediction
############################################
############################################
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K
y <- hybrid2$Yield
X1 <- model.matrix(~ Location, data = hybrid2);dim(X1)
Z1 <- model.matrix(~ GCA1 -1, data = hybrid2);dim(Z1)
Z2 <- model.matrix(~ GCA2 -1, data = hybrid2);dim(Z2)
Z3 <- model.matrix(~ SCA -1, data = hybrid2);dim(Z3)

K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)     
# Realized IBS relationships for set of parents 1
K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)     
# Realized IBS relationships for set of parents 2
S <- kronecker(K1, K2) ; dim(S)   
rownames(S) <- colnames(S) <- levels(hybrid2$SCA)
# Realized IBS relationships for cross 
#(as the Kronecker product of K1 and K2)

ETA <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=S))
# run the next line, ommited for CRAN time limitations
# ans <- mmer(y=y, X=X1, Z=ETA)
# ans$var.comp
# summary(ans)

#############################
## COMPARE WITH MCMCglmm ##
#############################
##== the same model run in MCMCglmm:
#library(MCMCglmm)
# pro <- list(GCA1 = as(solve(K1), "sparseMatrix"), GCA2 = as(solve(K2),
#      + "sparseMatrix"), SCA = as(solve(S), "sparseMatrix") )
#system.time(mox <- MCMCglmm(Yield ~ Location, random = ~ GCA1 + GCA2 + SCA,
#      + data = hybrid2, verbose = T, ginverse=pro))
## Takes 7:13 minutes in MCMCglmm, in sommer only takes 7 seconds


## it is also possible to do GWAS for hybrids, separatting and accounting
## for effects of GCA1, GCA2, SCA

#############################
## COMPARE WITH pedigreemm ##
#############################
# library(pedigreemm)
#A <- as.matrix(getA(pedCowsR))
#y <- milk$milk
#Z1 <- model.matrix(~id-1, data=milk); dim(Z1)
#vv <- match(unique(milk$id), gsub("id","",colnames(Z1)))
#K1<- A[vv,vv]; dim(K1) 
#Z2 <- model.matrix(~as.factor(herd)-1, data=milk); dim(Z2)
#ETA<- list(list(Z=Z1, K=K1),list(Z=Z2))
#fm3 <- mmer(y=y, Z=ETA) # or using mmer2 would look:
#fm3 <- mmer2(fixed=milk ~ 1, random = ~ id + herd, 
#             + G=list(id=K1), data=milk, draw=FALSE)
#summary(fm3)
# Try pedigreemm but takes longer, is an extension of lme4
#fm2 <- pedigreemm(milk ~ (1 | id) + (1 | herd),data = milk, pedigree = list(id= pedCowsR))
#plot(fm3$u.hat[[1]], ranef(fm2)$id[,1])
#plot(fm3$u.hat[[2]], ranef(fm2)$herd[,1])
# a big data frame with 3397 rows and 1359 animals analyzed
# pedigreemm takes 4 min, sommer takes 1 minute

#####################################
## PREDICTING SPECIFIC PERFORMANCE ##
## within biparental population    ##
#####################################

#data(CPdata)
#CPpheno <- CPdata$pheno
#CPgeno <- CPdata$geno
## look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
## fit a model including additive and dominance effects
#y <- CPpheno$color
#Za <- diag(length(y))
#Zd <- diag(length(y))
#A <- A.mat(CPgeno)
#D <- D.mat(CPgeno)

#y.trn <- y # for prediction accuracy
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#y.trn[ww] <- NA

### ADDITIVE MODEL ###
#ETA.A <- list(list(Z=Za,K=A))
#ans.A <- mmer(y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
### ADDITIVE-DOMINANT MODEL ###
#ETA.AD <- list(list(Z=Za,K=A),list(Z=Zd,K=D))
#ans.AD <- mmer(y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")

### 0.63 accuracy !!!! 4 percent increment!!

}