\name{sommer-package}
\alias{sommer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{Solving Mixed Model Equations in R}
\description{
The sommer package has been developed to provide R users with free code to understand how most common algorithms in mixed model analysis work related to genetic analysis and other general experiments, but at the same time allowing to perform their real analysis in diploid and polyploid organisms. This package allows the user to estimate variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects and obtain other parameters such as BLUPs, BLUEs, residuals, fitted values, variances for fixed and random effects, etc. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominant (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices that have been shown to increase prediction accuracy. The package provides flexibility to fit other genetic models such as full and half diallel models as well. 

Several datasets for learning how to use the sommer package are available; the \code{\link{HDdata}} dataset will help users know how to fit half and full diallel designs. The \code{\link{cornHybrid}} dataset contains data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset contains data to teach how to fit genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominant and epistatic effects. The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions. 

The core of the package is the function \code{\link{mmer}}, and the more user-friendly, but less flexible \code{\link{mmer2}}, both focused in solving mixed models of the form:

.

.............................    y = Xb + Zu + e  

.

b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, K*sigma(u)]  .......where: K*sigma(u) = G

e ~ N[0, I*sigma(e)]  .......where: I*sigma(e) = R

y ~ N[Xb, var(Zu+e)]   ......where;

var(y) = var(Zu+e) = ZGZ+R = V which is the phenotypic variance

.

The function allows the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with the following design:

.

fixed = only intercept.....................b ~ N[b.hat, 0]

random = GCA1 + GCA2 + SCA.................u ~ N[0, G]       

.

where G can have the form:

.

|K*sigma2(gca1).....................0..........................0.........|
    
|.............0.............S*sigma2(gca2).....................0.........| = G

|.............0....................0......................W*sigma2(sca)..|

.

And can be extended to any number of covariance structures. This would be specified in the \code{\link{mmer}} function as:

.

X1 <- matrix(1,length(y),1)      incidence matrix for intercept only

ETA1 <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=K3))      for 3 random effects

.

where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively and K1, K2, K3 their var-cov matrices. Now the fitted model will be:

.

ans <- mmer(y=y, X=X1, Z=ETA1)

.

FOR DETAILS ON HOW THE "AI", EM" AND "EMMA" ALGORITHMS WORK PLEASE REFER TO \code{\link{AI}} , \code{\link{EM}} AND \code{\link{EMMA}} help pages.

In addition, the package contains a very nice function to plot genetic maps with numeric variable or traits next to the LGs, see the \code{\link{map.plot2}} function to see how easy can be done. The package contains other functions:

\code{\link{transp}} function transform a vector of colors in transparent colors.

\code{\link{fdr}} calculates the false discovery rate for a vector of p-values.

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the \code{\link{mmer}} function which is the core of the sommer package.

}
\references{

Covarrubias-Pazaran G (2016) sommer: An R package for mixed model analysis education. R package version 1.1. URL https://cran.r-project.org/web/packages/sommer/.

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Henderson C.R. 1975. Best Linear Unbiased Estimation and Prediction under a Selection Model. Biometrics vol. 31(2):423-447.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

}
\keyword{ R package }
\seealso{
http://cggl.horticulture.wisc.edu/home-page/
}
\examples{

############################################
############################################
# breeding values with 1 variance component
############################################
############################################
#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- ifelse(runif(1000)<0.5,-1,1)
}
#random phenotypes
u <- rnorm(1000)
g <- as.vector(crossprod(t(M),u))
h2 <- 0.5  #heritability
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

Z1 <- diag(length(y))
ETA <- list( list(Z=Z1, K=A.mat(M)))
ans <- mmer(y=y, Z=ETA, method="EMMA") 
# change method to "AI" which is the default and more efficient method

############################################
############################################
# GWAS with 1 variance component and one A matrix
############################################
############################################
ETA <- list( list(Z=Z1, K=A.mat(M))) 
# random effects for genotypes, RUN IT:
# ans <- mmer(y=y, Z=ETA, W=M, method="EMMA") 

###################################################
###################################################
# IMPORT DATA FOR  ESTIMATING 3 VARIANCE COMPONENTS
###################################################
###################################################
# Import phenotypic data on inbred performance
## Full data
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K # genomic relationship matrix
############################################
############################################
# breeding values with 3 variance components
# hybrid prediction
############################################
############################################
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
# Realized IBS relationships for cross 
#(as the Kronecker product of K1 and K2)
rownames(S) <- colnames(S) <- levels(hybrid2$SCA)

ETA <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=S))
## run the next line, it was ommited for CRAN time limitations
#ans <- mmer(y=y, Z=ETA) 
#summary(ans)

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