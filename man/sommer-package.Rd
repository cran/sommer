\name{sommer-package}
\alias{sommer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{Solving Mixed Model Equations in R}
\description{


Multivariate-Univariate linear mixed model solver for multiple random effects allowing the specification of variance covariance structures. ML/REML estimates can be obtained using one of the following methods; the \strong{Direct-Inversion Newton-Raphson}, \strong{Direct-Inversion Average Information}, the \strong{MME-based Expectation-Maximization}, \strong{MME-based Average Information}, or the \strong{Efficient Mixed Model Association} algorithms. Designed for genomic prediction and genome wide association studies (GWAS) to include additive, dominance and epistatic relationship structures or other covariance structures, but also functional as a regular mixed model program. Multivariate models (multiple responses) can be fitted currently with any of the previously mentioned direct-inversion algorithms. In addition, \strong{rrBLUP results can be recreated selecting the EMMA method}.  

The sommer package has been developed to provide R users with open-source code to understand how most popular likelihood algorithms in mixed model analysis work, related to genetic analysis and other general experiments, but at the same time allowing to perform their real analysis in diploid and polyploid organisms with small and big data sets. The core of the package is the function \code{\link{mmer}} focused in solving multivariate linear mixed models by likelihood methods. \strong{There's 2 ways of fitting a model, specifying all matrices using \code{\link{mmer}} or using the formula-based solver \code{\link{mmer2}}}.This package allows the user to estimate variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects and obtain other parameters such as BLUPs, BLUEs, residuals, fitted values, variances-covariances for fixed and random effects, etc. The package is focused on genomic prediction (and genomic selection) and GWAS analysis, although general mixed models can be fitted as well. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices for diploid and polyploid organisms that have been shown to increase prediction accuracy under certain scenarios. The package provides flexibility to fit other genetic models such as full and half diallel models as well. In addition the \code{\link{pin}} function can be used to estimate standard errors for linear combinations of variance components (i.e. ratios).

The package has been equiped with several datasets to learn how to use the sommer package; the \code{\link{HDdata}} and \code{\link{FDdata}} datasets will introduce users to fit half and full diallel designs respectively. The \code{\link{h2}} dataset shows how to calculate heritability. The \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets contain data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset teaches how to do genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. The same data example is used to show how to include the top GWAS hits as fixed effects in the GBLUP model to increase prediction accuracy, the examples can be found in the \code{\link{hits}} documentation.  The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. The second example in \code{\link{Technow_data}} data shows how to perform GWAS in single cross hybrids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. The \code{\link{RICE}} dataset can teach users how to select the best training population using the \code{\link{TP.prep}} function in an applied breeding program, and we show comparisons with other methods. Traces of selection can be obtained using markers with the \code{\link{eigenGWAS}} function. A comparison of genomic prediction using univariate versus multivariate models is shown in the example#2 in the \code{\link{CPdata}} help page. Examples of multivariate models can be found in the example #3 of the \code{\link{CPdata}}. Furthermore, since v2.3 we have added the \code{\link{nna}} function which does nearest neighbor to create a new data frame adjusted for spatial variation. The \code{\link{ExpDesigns}} data contains several datasets to analyze experimental designs relevant to plant breeding and several detailed examples are available. Useful functions for analyzing such designs are included such as the \code{\link{blocker}} and \code{\link{fill.design}}. The \code{\link{gryphondata}} data contains an example of an animal model including pedigree information. The \code{\link{BTdata}} dataset contains an animal (birds) model to show the use of the \code{\link{pin}} function. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions.

\strong{For a short tutorial} of how to perform different genetic analysis in sommer \strong{please type}:

\strong{vignette("sommer")}

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the \code{\link{mmer}} function which is the core of the sommer package.  

Additional functions for genetic analysis have been included such as False Discovery Rate calculation (\code{\link{fdr}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}), phasing F1 or CP populations (\code{\link{phase.F1}}). We have also included the famous \code{\link{transp}} function to get transparent colors. 

\strong{Please update 'sommer' every month using:}


install.packages("sommer")


\strong{If you have good knowledge of residual structures (REML estimates of residual structures) and are willing to help me to implement them in sommer please send me an email.}

=======================================================
=======================================================

\strong{MODELS ENABLED}

The mmer function fits linear mixed models by likelihood methods allowing the used of known covariance structures for random effects. If not provided independence is conferred. This program focuses in the mixed model of the form:

'

\deqn{Y = X \beta + Z u + \epsilon}

'

with distributions:

'

\deqn{Y ~ MVN ( X\beta+Zu, var(Z u + \epsilon) )}

where;

'

\deqn{\beta  ~ N (\beta, 0)} 

\deqn{u  ~ N (0, G)} 

where G is equal to:

'

\tabular{rrrrr}{
  K1*\eqn{\sigma2}(u1)  \tab \tab 0 \tab \tab 0\cr
  0 \tab \tab K2*\eqn{\sigma2}(u2) \tab \tab 0\cr
  ... \tab \tab ... \tab \tab ...\cr
  0 \tab \tab 0 \tab \tab Ki*\eqn{\sigma2}(ui)\cr
}


'

for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and

'

\deqn{\epsilon  ~ N (0, R)}

where: \eqn{R = I * \sigma2 \epsilon}

'


and Var(Y) = Var(Zu+e) = ZGZ+R = V which is the phenotypic variance. Estimation of variance components are optimized with any of the 4 methods available; NR, AI, EM, EMMA (rrBLUP method).

=======================================================
=======================================================

\strong{GWAS MODELS}

The GWAS models in the sommer package are enabled by using the W argument, which is expected to be a numeric marker matrix. Markers are treated as fixed effects according to the model proposed by Yu et al. (2006) for diploids, and Rosyara et al. (2016) (for polyploids). The matrices X and W are both fixed effects, but they are separated by 2 different arguments to distinguish factors such as environmental and design factors for the argument "X" and markers with "W".
  
  The genome-wide association analysis is based on the mixed model:

\deqn{y = X \beta + Z g + W \tau + e}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The argument "P3D" was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

Multivariate GWAS are based in Covarrubias-Pazaran et al. (2017, Submitted), which adjusts betas for all response variables and then does the regular GWAS with such adjusted betas or marker effects.

For extra details about the methods please read the canonical papers listed in the References section.

=======================================================
=======================================================

As an example, imagine a mixed model with three random effects G takes the form:

'

\tabular{rrrrr}{
  K1*\eqn{\sigma2}(gca1)  \tab 0 \tab 0\cr
  0 \tab K2*\eqn{\sigma2}(gca2) \tab 0 \tab  =  G\cr
  0 \tab 0 \tab K3*\eqn{\sigma2}(sca)\cr
}

'

This mixed model would be specified in the \code{\link{mmer}} function as:

'

X1 <- matrix(1,length(y),1)     # incidence matrix for intercept only

ETA <- list(

gca1=list(Z=Z1, K=K1), 

gca2=list(Z=Z2, K=K2), 

sca=list(Z=Z3, K=K3)

)      #for 3 random effects

'

where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively created using the \code{\link{model.matrix}} function and K1, K2, K3 are their var-cov matrices. Now the fitted model will be typed as:

'

ans <- mmer(Y=Y, X=X1, Z=ETA)

or


'

or if a data frame is available:

'

ans <- mmer2(Y~1, random= ~ g(gca1) + g(gca2) + g(sca), 
             G=list(gca1=K1, gca2=K2, sca=K3), 
             data=yourdata)



=======================================================
=======================================================

The mmer function has been enhanced by adding the argument EIGEND which allows an eigen decomposition of the additive relationship matrix to accelerate genomic prediction models in the order of hundreds of times compared to classical EMMA or AI when dense covariance structures.

Finally, feel free to get in touch with me if you have any questions or suggestions at: 

cova_ruber@live.com.mx


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
https://scholar.google.com/citations?user=Ic0qn-kAAAAJ&hl=en
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
colnames(M) <- paste("geno",1:dim(M)[2], sep="")
rownames(M) <- paste("mark",1:dim(M)[1], sep="")
####=========================================####
#### simulate phenotypes
####=========================================####
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
M <- t(M); M[1:5,1:5]

####=========================================####
#### fit the model using mmer (matrix based)
####=========================================####
Z1 <- diag(length(y))
K1 <-A.mat(M)
ETA <- list( add=list(Z=Z1, K=K1) )
#ans <- mmer(Y=y, Z=ETA)
#summary(ans)

#### run the same but as GWAS 
#### just add the marker matrix in the argument W
#### markers are fixed effects

#ans <- mmer(Y=y, Z=ETA, W=M)
#summary(ans)

####=========================================####
#### fit the model using mmer2 (formula-based)
####=========================================####
K<-A.mat(M) # additive relationship matrix
dat <- data.frame(y=y,id=rownames(M));head(dat)
#ans <- mmer2(y~1, random=~g(id), G=list(id=K), data=dat)
#summary(ans)

#### run the same but as GWAS 
#### just add the marker matrix in the argument W
#### markers are fixed effects

#ans <- mmer2(y~1, random=~g(id), G=list(id=K), W=M, data=dat)
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
####=========================================####
####=========================================####


####=========================================####
#### Hybrid prediction using mmer (matrix-based)
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


####=========================================####
#### Hybrid prediction using mmer2 (formula-based)
####=========================================####
# data(cornHybrid)
# hybrid2 <- cornHybrid$hybrid # extract cross data
# A <- cornHybrid$K
# ####=========================================####
# #### Realized IBS relationships for set of parents 1
# ####=========================================####
# K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1) 
# ####=========================================####
# #### Realized IBS relationships for set of parents 2
# ####=========================================####
# K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)
# ####=========================================####
# #### Realized IBS relationships for cross 
# #### (as the Kronecker product of K1 and K2)
# ####=========================================####
# S <- kronecker(K1, K2) ; dim(S)   
# rownames(S) <- colnames(S) <- levels(hybrid2$SCA)
# 
# ans <- mmer2(Yield~Location, random=~g(GCA1)+ g(GCA2) + g(SCA),
#              G=list(GCA1=K1, GCA2=K2, SCA=S), data=hybrid2)
# summary(ans)

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

####=========================================####
#### using mmer (matrix-based)
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
#A <- A.mat(CPgeno) # additive relationship matrix
#D <- D.mat(CPgeno) # dominant relationship matrix

####================####
#### ADDITIVE MODEL ####
####================####
#ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- mmer(Y=y, Z=ETA.A)
#summary(ans.A)
####=========================####
#### ADDITIVE-DOMINANCE MODEL ####
####=========================####
#ETA.AD <- list(add=list(Z=Za,K=A), dom=list(Z=Zd,K=D))
#ans.AD <- mmer(Y=y, Z=ETA.AD)
#summary(ans.AD)
### greater accuracy !!!! 4 percent increment!!
### we run 100 iterations, 4 percent increment in general

####=========================================####
#### same using mmer2 (formula-based)
####=========================================####

# data(CPdata)
# CPpheno <- CPdata$pheno
# CPgeno <- CPdata$geno
# ### look at the data
# head(CPpheno); CPgeno[1:5,1:5]
# ###=========================================####
# ### fit a model including additive and dominance effects
# ###=========================================####
# A <- A.mat(CPgeno)
# D <- D.mat(CPgeno)
# ####=========================================####
# #### ADDITIVE MODEL 
# ####=========================================####
# ans.A <- mmer2(Yield~1,random=~ g(id), G=list(id=A), 
#                data=CPpheno)
# summary(ans.A)
# ####=========================================####
# #### ADDITIVE-DOMINANCE MODEL 
# ####=========================================####
# CPpheno$idd <- CPpheno$id 
# ans.AD <- mmer2(Yield~1,random=~ g(id) + g(idd), 
#                 G=list(id=A, idd=D), data=CPpheno)
# summary(ans.AD)

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

# data(CPdata)
# CPpheno <- CPdata$pheno; CPpheno$idd <- CPpheno$id 
# CPgeno <- CPdata$geno
# ###=========================================####
# ### Do incidence and relationship matrices
# ###=========================================####
# A <- A.mat(CPgeno)
# D <- D.mat(CPgeno)
# ###=========================================####
# ### Only additive model using all traits
# ###=========================================####
# head(CPpheno)
# ans.A <- mmer2(cbind(color,Yield,FruitAver)~1, random = ~ g(id),
#               G=list(id=A), n.cores = 3, data=CPpheno)
# summary(ans.A)
# ###=========================================####
# ### Additive + dominance model
# ###=========================================####
# head(CPpheno)
# ans.AD <- mmer2(cbind(color,Yield,FruitAver)~1, random=~g(id)+g(idd),
#                G=list(id=A, idd=D), n.cores = 3, data=CPpheno)
# summary(ans.AD)

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
# data(CPdata)
# CPpheno <- CPdata$pheno
# CPgeno <- CPdata$geno
# ### look at the data
# head(CPpheno)
# CPgeno[1:5,1:5]
# ## fit a model including additive effects
# A <- A.mat(CPgeno) 
# ### MAKE SURE YOU SET 'MVM=TRUE'
# ans.A <- mmer2(cbind(color,Yield,FruitAver)~1, random = ~ g(id),
#               G=list(id=A),MVM=TRUE, data=CPpheno)
# summary(ans.A)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 6
#### OVERLAY IN MODELS
####=========================================####
####=========================================####


####=========================================####
#### using mmer (matrix-based)
####=========================================####
data(HDdata)
head(HDdata)

#### GCA matrix for half diallel using male and female columns

Z1 <- overlay(HDdata[,c("male","female")])

#### SCA matrix

Z2 <- model.matrix(~as.factor(geno)-1, data=HDdata)

#### Define response variable and run

y <- HDdata$sugar
ETA <- list(GCA=list(Z=Z1), SCA=list(Z=Z2)) # Zu component
modHD <- mmer(Y=y, Z=ETA)
summary(modHD)

####=========================================####
#### using overlay with mmer2 function (formula-based)
####=========================================####
# data(HDdata)
# head(HDdata)
# HDdata$female <- as.factor(HDdata$female)
# HDdata$male <- as.factor(HDdata$male)
# HDdata$geno <- as.factor(HDdata$geno)
# #### model using overlay
# modh <- mmer2(sugar~1, random=~female + and(male) + geno, 
#               data=HDdata)
# summary(modh)

# #### model using overlay [and(.)] and covariance structures [g(.)]

# A <- diag(7); A[1,2] <- 0.5; A[2,1] <- 0.5 # fake covariance structure
# colnames(A) <- as.character(1:7); rownames(A) <- colnames(A);A
# 
# modh2 <- mmer2(sugar~1, random=~g(female) + and(g(male)) + geno, 
#               G=list(female=A, male=A),data=HDdata)
# summary(modh2)

# data(HDdata)
# head(HDdata)
# HDdata$female <- as.factor(HDdata$female)
# HDdata$male <- as.factor(HDdata$male)
# HDdata$geno <- as.factor(HDdata$geno)
# #### model using overlay
# modh <- mmer2(sugar~1, random=~female + and(male) + geno, 
#               data=HDdata)
# summary(modh)

}