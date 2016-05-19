\name{sommer-package}
\alias{sommer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{Solving Mixed Model Equations in R}
\description{
The sommer package has been developed to provide R users with open-source code to understand how most common algorithms in mixed model analysis work related to genetic analysis and other general experiments, but at the same time allowing to perform their real analysis in diploid and polyploid organisms with small and big data sets. This package allows the user to estimate variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects and obtain other parameters such as BLUPs, BLUEs, residuals, fitted values, variances for fixed and random effects, etc. The package is focused on genomic prediction (or genomic selection) and GWAS analysis although general mixed models can be fitted as well. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices that have been shown to increase prediction accuracy. The package provides flexibility to fit other genetic models such as full and half diallel models as well. 

Several datasets for learning how to use the sommer package are available; the \code{\link{HDdata}} and \code{\link{FDdata}} datasets will help users know how to fit half and full diallel designs respectively. The \code{\link{h2}} dataset shows how to calculate heritability. The \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets contain data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset contains data to teach how to fit genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominant and epistatic effects. The same data example shows how to do what has been denominated bagging GBLUP which increases prediction accuracy and examples can be found in the \code{\link{bag}} documentation.  The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. The second example in \code{\link{Technow_data}} data shows how to perform GWAS in single cross hybrids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. The \code{\link{RICE}} dataset can teach users how to select the best training population using the \code{\link{TP.prep}} function in an applied breeding program, and we show comparison some comparison with other methods. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions. 

For a short tutorial of how to perform different genetic analysis in sommer please type:

vignette("sommer")

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the \code{\link{mmer}} function which is the core of the sommer package. The core algorithm in our mixed model solver is the Direct Average Information proposed by Lee et al. (2015) which surpasses in performance the MME-based Average Information (Gilmour et al. 1995) when dense covariance structures are present (GBLUP and GWAS case). In addition the Expectation-Maximization and Newton-Raphson algorithms are also available and can deal with multiple variance components with var-cov structures. When the variance covariance matrices are identity matrices (sparse), or longitudinal data is being analyzed we recommend shifting to use the lme4 package for such cases of NO covariance structures. 

The core of the package is the function \code{\link{mmer}}, and the more user-friendly (asreml sintax type) but less flexible \code{\link{mmer2}}, both focused in solving mixed models of the form:

.

.............................    y = Xb + Zu + e   ........ with distributions:

.

y ~ N[Xb, var(Zu+e)]   ......where;

.

b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, G]  ....... where G is equal to:

.

|K1*sigma2(u1)......................0...........................0.........|
    
|.............0.............K2*sigma2(u2).......................0.........| = G

|..........................................................................................|

|.............0....................0.........................Ki*sigma2(ui)...|

.

for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and

.

e ~ N[0, R]  .....................where: I*sigma(e) = R

.

also Var(y) = Var(Zu+e) = ZGZ+R = V which is the phenotypic variance

-------------------------------------------------------------------------------

The functions in the sommer packages allow the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with three random effects:

.

fixed = only intercept...................................b ~ N[b.hat, 0]

random = GCA1 + GCA2 + SCA.................u ~ N[0, G]       

.

then G takes the form:

.

|K1*sigma2(gca1).....................0..........................0.........|
    
|.............0.............K2*sigma2(gca2).....................0.........| = G

|.............0....................0......................K3*sigma2(sca)...|

.

This mixed model would be specified in the \code{\link{mmer}} function as:

.

X1 <- matrix(1,length(y),1)      incidence matrix for intercept only

ETA <- list(gca1=list(Z=Z1, K=K1), gca2=list(Z=Z2, K=K2), sca=list(Z=Z3, K=K3))      for 3 random effects

.

where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively created using the \code{\link{model.matrix}} function and K1, K2, K3 are their var-cov matrices. Now the fitted model will be typed as:

.

ans <- mmer(y=y, X=X1, Z=ETA)

or

ans <- mmer2(y~1, random= ~ gca1 + gca2 + sca, G=list(gca1=K1, gca2=K2, sca=K3), data=yourdata)

-------------------------------------------------------------------------------------

The package uses Average Information \code{\link{AI}} , Expectation Maximization \code{\link{EM}}, Newton-Raphson \code{\link{NR}}, or Efficient Mixed Model Association \code{\link{EMMA}} algorithms. I have included those functions as independent in case you want to check how they work internally.

In addition, the package contains a very nice function to plot genetic maps with numeric variable or traits next to the LGs, see the \code{\link{map.plot}} function to see how easy can be done. The package contains other functions:

\code{\link{transp}} function transform a vector of colors in transparent colors.

\code{\link{fdr}} calculates the false discovery rate for a vector of p-values.

The mmer function has been enhanced by adding the argument MTG2 which allows an eigen decomposition of the additive relationship matrix to accelerate genomic prediction models in the order of tens of thousands of individuals.

Finally, feel free to get in touch with me if you have any questions or suggestions at: 

covarrubiasp@wisc.edu

I'll be glad to help or answer any question. We have spent a valuable amount of time developing this package. Please cite us in your publication. Type 'citation("sommer")' to know how to cite it.


}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. https://cran.rstudio.com/web/packages/sommer/ 

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Henderson C.R. 1975. Best Linear Unbiased Estimation and Prediction under a Selection Model. Biometrics vol. 31(2):423-447.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

Abdollahi Arpanahi R, Morota G, Valente BD, Kranis A, Rosa GJM, Gianola D. 2015. Assessment of bagging GBLUP for whole genome prediction of broiler chicken traits. Journal of Animal Breeding and Genetics 132:218-228.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

}
\keyword{ R package }
\seealso{
http://cggl.horticulture.wisc.edu/home-page/
}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
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
ETA <- list( list(Z=Z1, K=A.mat(M)))
#ans <- mmer(y=y, Z=ETA, method="EMMA")
#summary(ans)

####=========================================####
#### run the same but as GWAS 
#### just add the marker matrix in the argument W
#### markers are fixed effects
####=========================================####

#ans <- mmer(y=y, Z=ETA, W=M, method="EMMA")
#summary(ans)

############################################
############################################
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

#ETA <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=S))
#ans <- mmer(y=y, X=X1, Z=ETA)
#ans$var.comp
#summary(ans)

############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 3
#### PREDICTING SPECIFIC PERFORMANCE
#### within biparental population
####=========================================####
####=========================================####

#data(CPdata)
#CPpheno <- CPdata$pheno
#CPgeno <- CPdata$geno
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
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for one fifth of the population
#y.trn[ww] <- NA

####=========================================####
#### ADDITIVE MODEL
####=========================================####
#ETA.A <- list(list(Z=Za,K=A))
#ans.A <- mmer(y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================================####
#### ADDITIVE-DOMINANCE MODEL 
####=========================================####
#ETA.AD <- list(list(Z=Za,K=A),list(Z=Zd,K=D))
#ans.AD <- mmer(y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
#summary(ans.A)
#summary(ans.AD)
### 0.63 accuracy !!!! 4 percent increment!!
}