\name{sommer-package}
\alias{sommer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{Solving Mixed Model Equations in R}
\description{

Mixed model equation solver allowing the specification of variance covariance structures for random effects. ML/REML estimates are obtained using the \strong{Average Information}, \strong{Expectation-Maximization}, \strong{Newton-Raphson}, or \strong{Efficient Mixed Model Association} algorithms. Designed for genomic prediction and genome wide association studies (GWAS) to include additive, dominance and epistatic relationship structures or other covariance structures in R, but also functional as a regular mixed model program. Multivariate models are enabled for a single random effect.

The sommer package has been developed to provide R users with open-source code to understand how most popular likelihood algorithms in mixed model analysis work, related to genetic analysis and other general experiments, but at the same time allowing to perform their real analysis in diploid and polyploid organisms with small and big data sets. The core of the package is the function \code{\link{mmer}}, and the more user-friendly, \code{\link{mmer2}}, both focused in solving linear mixed models by likelihood methods. This package allows the user to estimate variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects and obtain other parameters such as BLUPs, BLUEs, residuals, fitted values, variances for fixed and random effects, etc. The package is focused on genomic prediction (or genomic selection) and GWAS analysis, although general mixed models can be fitted as well. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices for diploid and polyploid organisms that have been shown to increase prediction accuracy under certain scenarios. The package provides flexibility to fit other genetic models such as full and half diallel models as well. 

The package has been equiped with several datasets to learn how to use the sommer package; the \code{\link{HDdata}} and \code{\link{FDdata}} datasets will introduce users to fit half and full diallel designs respectively. The \code{\link{h2}} dataset shows how to calculate heritability. The \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets contain data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset teaches how to do genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. The same data example is used to show how to include the top GWAS hits as fixed effects in the GBLUP model to increase prediction accuracy, the examples can be found in the \code{\link{hits}} documentation.  The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. The second example in \code{\link{Technow_data}} data shows how to perform GWAS in single cross hybrids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. The \code{\link{RICE}} dataset can teach users how to select the best training population using the \code{\link{TP.prep}} function in an applied breeding program, and we show comparison some comparison with other methods. Traces of selection can be obtained using markers with the \code{\link{eigenGWAS}} function. A comparison of genomic prediction using univariate versus multivariate models is shown in the example#2 in the \code{\link{CPdata}} help page. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions. 

For a short tutorial of how to perform different genetic analysis in sommer please type:

vignette("sommer")

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the \code{\link{mmer}} function which is the core of the sommer package. The core algorithm in our mixed model solver is the Direct Average Information proposed by Lee et al. (2015) which surpasses in performance the MME-based Average Information (Gilmour et al. 1995) when dense covariance structures are present (GBLUP and GWAS case). In addition the Expectation-Maximization and Newton-Raphson algorithms are also available and can deal with multiple variance components with var-cov structures, in addition to the popular EMMA algorithm. When the variance covariance matrices are identity matrices (sparse), or longitudinal data is being analyzed we recommend shifting to \strong{use the lme4 package for} such \strong{cases of NO covariance structures}.  

Additional functions for genetic analysis have been included such as False Discovery Rate calculation (\code{\link{fdr}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}). We have also included the famous \code{\link{transp}} function to get transparent colors. 

PLEASE UPDATE 'sommer' EVERY MONTH USING:

'
install.packages("sommer")

The core of the package is the function \code{\link{mmer}}, and the more user-friendly (asreml sintax type) \code{\link{mmer2}}, both focused in solving mixed models of the form:

\deqn{Y = X \beta + Z u + \varepsilon}  with distributions:


\deqn{Y - MVN(X \beta, Var(Z u + \varepsilon))}


where;


\deqn{u - MVN(0, G)}


where G is equal to:


\deqn{K_1 \sigma^2_1 --- 0 --- 0}
\deqn{0 ---- K_2 \sigma^2_2 -- 0}
\deqn{...................................}
\deqn{0 ---- 0 -- K_i \sigma^2_i}


for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and


\deqn{\epsilon - MVN(0, R)}  

where:

\deqn{R = I \sigma^2_e}

also

\deqn{Var(y) = Var(Zu + \epsilon) = ZGZ' + R = V}  

which is the phenotypic variance


'

The functions in the sommer packages allow the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with three random effects:



fixed = only intercept  

\deqn{\beta - MVN(\beta, 0)}


random = GCA1 + GCA2 + SCA  

\deqn{u - MVN(0, G)}       



then G takes the form:



\deqn{K_1 \sigma^2_1 --- 0 --- 0}
\deqn{0 ---- K_2 \sigma^2_2 -- 0}
\deqn{0 ----- 0 -- K_3 \sigma^2_3}


This mixed model would be specified in the \code{\link{mmer}} function as:

# incidence matrix for intercept only


X1 <- matrix(1,length(y),1)     


# for 3 random effects


ETA <- list(

 gca1=list(Z=Z1, K=K1),

           gca2=list(Z=Z2, K=K2),

           sca=list(Z=Z3, K=K3)

)      


where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively created using the \code{\link{model.matrix}} function and K1, K2, K3 are their var-cov matrices. Now the fitted model will be typed as:

'

ans <- mmer(Y=y, X=X1, Z=ETA)

'

or if a data frame is available:

'

ans <- mmer2(Y~1, random= ~ gca1 + gca2 + sca, 
             G=list(gca1=K1, gca2=K2, sca=K3), 
             data=yourdata)


'


The package uses Average Information \code{\link{AI}} , Expectation Maximization \code{\link{EM}}, Newton-Raphson \code{\link{NR}}, or Efficient Mixed Model Association \code{\link{EMMA}} algorithms. I have included those functions as independent in case you want to check how they work internally.

The mmer function has been enhanced by adding the argument EIGEND which allows an eigen decomposition of the additive relationship matrix to accelerate genomic prediction models in the order of hundreds of times compared to classical EMMA or AI.

Finally, feel free to get in touch with me if you have any questions or suggestions at: 

covarrubiasp@wisc.edu


I'll be glad to help or answer any question. We have spent a valuable amount of time developing this package. Please cite us in your publication. Type 'citation("sommer")' to know how to cite it.


}
\references{

Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

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
ETA <- list(add=list(Z=Z1, K=A.mat(M)))
#ans <- mmer(Y=y, Z=ETA, method="EMMA")
#summary(ans)

####=========================================####
#### run the same but as GWAS 
#### just add the marker matrix in the argument W
#### markers are fixed effects
####=========================================####

#ans <- mmer(Y=y, Z=ETA, W=M, method="EMMA")
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

#ETA <- list(GCA1=list(Z=Z1, K=K1), GCA2=list(Z=Z2, K=K2), SCA=list(Z=Z3, K=S))
#ans <- mmer(Y=y, X=X1, Z=ETA)
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
#ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- mmer(Y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================================####
#### ADDITIVE-DOMINANCE MODEL 
####=========================================####
#ETA.AD <- list(add=list(Z=Za,K=A),dom=list(Z=Zd,K=D))
#ans.AD <- mmer(Y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
#summary(ans.A)
#summary(ans.AD)
### 0.63 accuracy !!!! 4 percent increment!!


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
#### EXAMPLE 4
#### PARALLEL RESPONSES
#### accelerating analysis
#### using 'n.cores' argument for parallelization
####=========================================####
####=========================================####

#data(CPdata)
#CPpheno <- CPdata$pheno
#CPgeno <- CPdata$geno
####=========================================####
#### Calculate incidence and relationship matrices
####=========================================####
#N <- dim(CPpheno)[1]
#Za <- diag(N)
#Zd <- diag(N)
#A <- A.mat(CPgeno)
#D <- D.mat(CPgeno)
####=========================================####
#### Only additive model using all traits
####=========================================####
#head(CPpheno)
#ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- mmer(Y=CPpheno, Z=ETA.A, method="EMMA", n.cores = 4)
#summary(ans.A)
####=========================================####
#### Additive + dominance model
####=========================================####
#head(CPpheno)
#ETA.AD <- list(add=list(Z=Za,K=A),dom=list(Z=Zd,K=D))
#ans.AD <- mmer(Y=CPpheno, Z=ETA.AD, method="AI", n.cores = 4)
#summary(ans.AD)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 5
#### MULTIVARIATE MODEL
####=========================================####
####=========================================####
data(CPdata)
CPpheno <- CPdata$pheno
CPgeno <- CPdata$geno
### look at the data
head(CPpheno)
CPgeno[1:5,1:5]
## fit a model including additive and dominance effects
Y <- CPpheno
Za <- diag(dim(Y)[1])
A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
ETA.A <- list(add=list(Z=Za,K=A))
#### MAKE SURE YOU SET 'MVM=TRUE'
#ans.A <- mmer(Y=Y, Z=ETA.A, MVM=TRUE)
#summary(ans.A)

}