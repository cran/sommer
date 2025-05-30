\name{GWAS}
\alias{GWAS}
\title{

Genome wide association study analysis

}
\description{

THIS FUNCTION IS DEPRECATED. Fits a multivariate/univariate linear mixed model GWAS by likelihood methods (REML), see the Details section below. It uses the \code{\link{mmer}} function and its core coded in C++ using the Armadillo library to optimize dense matrix operations common in the derect-inversion algorithms. After the model fit extracts the inverse of the phenotypic variance matrix to perform the association test for the "p" markers. Please check the Details section (Model enabled) if you have any issue with making the function run.

The sommer package is updated on CRAN every 3-months due to CRAN policies but you can find the latest source at https://github.com/covaruber/sommer . This can be easily installed typing the following in the R console:

library(devtools)

install_github("covaruber/sommer")

This is recommended since bugs fixes will be immediately available in the GitHub source. \strong{For tutorials} on how to perform different analysis with sommer please look at the vignettes by typing in the terminal:

\strong{vignette("v1.sommer.quick.start")} 

\strong{vignette("v2.sommer.changes.and.faqs")}

\strong{vignette("v3.sommer.qg")}

\strong{vignette("v4.sommer.gxe")}

or visit \strong{https://covaruber.github.io}

}
\usage{

GWAS(fixed, random, rcov, data, weights, W,
    nIters=20, tolParConvLL = 1e-03, tolParInv = 1e-06, 
    init=NULL, constraints=NULL,method="NR", 
    getPEV=TRUE,naMethodX="exclude",
    naMethodY="exclude",returnParam=FALSE, 
    dateWarning=TRUE,date.warning=TRUE,verbose=FALSE,
    stepWeight=NULL, emWeight=NULL,
    M=NULL, gTerm=NULL, n.PC = 0, min.MAF = 0.05, 
    P3D = TRUE)
}
\arguments{
  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:
  
  \emph{response ~ covariate} for univariate models
  
  \emph{cbind(response.i,response.j) ~ covariate} for multivariate models
  
  The fcm function can be used to constrain fixed effects in multi-response models.
  
  }
  
  \item{random}{a formula specifying the name of the \strong{random effects}, i.e. \emph{random= ~ genotype + year}.
  
  Useful functions can be used to fit heterogeneous variances and other special models (\emph{see 'Special Functions' in the Details section for more information}):
  
  \code{\link{vsr}(...,Gu,Gt,Gtc)} is the main function to specify variance models and special structures for random effects. On the ... argument you provide the unknown variance-covariance structures (i.e. usr,dsr,at,csr) and the random effect where such covariance structure will be used (the random effect of interest). Gu is used to provide known covariance matrices among the levels of the random effect, Gt initial values and Gtc for constraints. Auxiliar functions for building the variance models are:
  
  }
  \item{rcov}{a formula specifying the name of the \strong{error term}, i.e. \emph{rcov= ~ units}.
  
  The functions that can be used to fit heterogeneous residual variances are the same used on the random term but the random effect is always "units", i.e. \emph{rcov=~vsr(dsr(Location),units)}
  
  }
  \item{data}{a data frame containing the variables specified in the formulas for response, fixed, and random effects.}
  \item{weights}{name of the covariate for weights. To be used for the product R = Wsi*R*Wsi, where * is the matrix product, Wsi is the square root of the inverse of W and R is the residual matrix.}
  \item{W}{Alternatively, instead of providing a vector of weights the user can specify an entire W matrix (e.g., when covariances exist). To be used first to produce Wis = solve(chol(W)), and then calculate R = Wsi*R*Wsi.t(), where * is the matrix product, and R is the residual matrix. Only one of the arguments weights or W should be used. If both are indicated W will be given the preference.}
  \item{nIters}{Maximum number of iterations allowed.  Default value is 15.}
  \item{tolParConvLL}{Convergence criteria.}
  \item{tolParInv}{tolerance parameter for matrix inverse used when singularities are encountered.}
  \item{init}{initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values for ALL var-cov components this argument is functional. It has to be provided as a list or an array, where each list element is one variance component and if multitrait model is pursued each element of the list is a matrix of variance covariance components among traits. Initial values can also be provided in the Gt argument of the \link{vsr} function.Is highly encouraged to use the Gt and Gtc arguments of the \link{vsr} function instead of this argument}
  \item{constraints}{when initial values are provided these have to be accompanied by their constraints. See the \link{vsr} function for more details on the constraints. Is highly encouraged to use the Gt and Gtc arguments of the \link{vsr} function instead of this argument.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. Direct-inversion Newton-Raphson \strong{NR} and Average Information \strong{AI} (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2015).}
  \item{getPEV}{a TRUE/FALSE value indicating if the program should return the predicted error variance and variance for random effects. This option is provided since this can take a long time for certain models where p > n by a big extent.}
   \item{naMethodX}{one of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the X matrices for fixed effects with the median value. If "exclude" is selected it will get rid of all rows with missing values for the X (fixed) covariates. The default is "exclude". The "include" option should be used carefully.}
   \item{naMethodY}{one of the three possible values; "include", "include2" or "exclude". If "include" is selected then the function will impute the response variables with the median value. The difference between "include" and "include2" is only available in the multitrait models when the imputation can happen for the entire matrix of responses or only for complete cases ("include2"). If "exclude" is selected it will get rid of rows in responses where missing values are present for the estimation of variance components. The default is "exclude".}
   \item{returnParam}{a TRUE/FALSE value to indicate if the program should return the parameters used for modeling without fitting the model.}
   \item{dateWarning}{a TRUE/FALSE value to indicate if the program should warn you when is time to update the sommer package.}
  \item{date.warning}{a TRUE/FALSE value to indicate if the program should warn you when is time to update the sommer package.}
  \item{verbose}{a TRUE/FALSE value to indicate if the program should return the progress of the iterative algorithm.}
    \item{stepWeight}{A vector of values (of length equal to the number of iterations) indicating the weight used to multiply the update (delta) for variance components at each iteration. If NULL the 1st iteration will be multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument can help to avoid that variance components go outside the parameter space in the initial iterations which doesn't happen very often with the NR method but it can be detected by looking at the behavior of the likelihood. In that case you may want to give a smaller weight to the initial 8-10 iterations.}
  \item{emWeight}{A vector of values (of length equal to the number of iterations) indicating with values between 0 and 1 the weight assigned to the EM information matrix. And the values 1 - emWeight will be applied to the NR/AI information matrix to produce a joint information matrix. If NULL weights for EM information matrix are zero and 1 for the NR/AI information matrix.}
  \item{M}{The marker matrix containing the marker scores for each level of the random effect selected in the gTerm argument, coded as numeric based on the number of reference alleles in the genotype call, e.g. (-1,0,1) = (aa,Aa,AA), levels in diploid individuals. Individuals in rows and markers in columns. No additional columns should be provided, is a purely numerical matrix. Similar logic applies to polyploid individuals, e.g. (-3,-2,-1,0,1,2,3) = (aaaa,aaaA,aaAA,Aaaa,AAaa,AAAa,AAAA). }
  \item{gTerm}{a character vector indicating the random effect linked to the marker matrix M (i.e. the genetic term) in the model. The random effect selected should have the same number of levels than the number of rows of M. When fitting only a random effect without a special covariance structure (e.g., dsr, usr, etc.) you will need to add the call 'u:' to the name of the random effect given the behavior of the naming rules of the solver when having a simple random effect without covariance structure. }
  \item{n.PC}{Number of principal components to include as fixed effects. Default is 0 (equals K model).}
  \item{min.MAF}{Specifies the minimum minor allele frequency (MAF). If a marker has a MAF less than min.MAF, it is assigned a zero score.}
  \item{P3D}{When P3D=TRUE, variance components are estimated by REML only once, without any markers in the model and then a for loop for hypothesis testing is performed. When P3D=FALSE, variance components are estimated by REML for each marker separately. The latter can be quite time consuming. As many models will be run as number of marker.}
    
}

\details{

\strong{Citation}

Type \emph{citation("sommer")} to know how to cite the sommer package in your publications. 

\strong{Models Enabled}

For details about the models enabled and more information about the covariance structures please check the help page of the package (\code{\link{sommer}}). In general the GWAS model implemented in sommer to obtain marker effect is a generalized linear model of the form:

b = (X'V-X)X'V-y

with X = ZMi

where:
b is the marker effect (dimensions 1 x mt)
y is the response variable (univariate or multivariate) (dimensions 1 x nt)
V- is the inverse of the phenotypic variance matrix (dimensions nt x nt)
Z is the incidence matrix for the random effect selected (gTerm argument) to perform the GWAS (dimensions nt x ut)
Mi is the ith column of the marker matrix (M argument) (dimensions u x m)

for t traits, n observations, m markers and u levels of the random effect. Depending if P3D is TRUE or FALSE the V- matrix will be calculated once and used for all marker tests (P3D=TRUE) or estimated through REML for each marker (P3D=FALSE).

 
\emph{vignette('sommer.start')}

\strong{Bug report and contact}

If you have any technical questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com. 

If you have any bug report please go to https://github.com/covaruber/sommer or send me an email to address it asap.

}

\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{Vi}{the inverse of the phenotypic variance matrix V^- = (ZGZ+R)^-1}
\item{sigma}{a list with the values of the variance-covariance components with one list element for each random effect.}
\item{sigma_scaled}{a list with the values of the scaled variance-covariance components with one list element for each random effect.}
\item{sigmaSE}{Hessian matrix containing the variance-covariance for the variance components. SE's can be obtained taking the square root of the diagonal values of the Hessian.}
\item{Beta}{a data frame for trait BLUEs (fixed effects).}
\item{VarBeta}{a variance-covariance matrix for trait BLUEs}
\item{U}{a list (one element for each random effect) with a data frame for trait BLUPs.}
\item{VarU}{a list (one element for each random effect) with the variance-covariance matrix for trait BLUPs.}
\item{PevU}{a list (one element for each random effect) with the predicted error variance matrix for trait BLUPs.}
\item{fitted}{ Fitted values y.hat=XB}
\item{residuals}{Residual values e = Y - XB}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{scores}{ marker scores (-log_(10)p) for the traits}
\item{method}{The method for extimation of variance components specified by the user.}
\item{constraints}{contraints used in the mixed models for the random effects.}

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   

Covarrubias-Pazaran G. 2018. Software update: Moving the R package sommer to multivariate mixed models for genome-assisted prediction. doi: https://doi.org/10.1101/354639

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

Rodriguez-Alvarez, Maria Xose, et al. Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics 23 (2018): 52-71.

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
#####========================================####
##### potato example
#####========================================####
# 
# data(DT_polyploid)
# DT <- DT_polyploid
# GT <- GT_polyploid
# MP <- MP_polyploid
# ####=========================================####
# ####### convert markers to numeric format
# ####=========================================####
# numo <- atcg1234(data=GT, ploidy=4);
# numo$M[1:5,1:5];
# numo$ref.allele[,1:5]
# 
# ###=========================================####
# ###### plants with both genotypes and phenotypes
# ###=========================================####
# common <- intersect(DT$Name,rownames(numo$M))
# 
# ###=========================================####
# ### get the markers and phenotypes for such inds
# ###=========================================####
# marks <- numo$M[common,]; marks[1:5,1:5]
# DT2 <- DT[match(common,DT$Name),];
# DT2 <- as.data.frame(DT2)
# DT2[1:5,]
# 
# ###=========================================####
# ###### Additive relationship matrix, specify ploidy
# ###=========================================####
# A <- A.mat(marks)
# ###=========================================####
# ### run it as GWAS model
# ###=========================================####
# ans2 <- GWAS(tuber_shape~1,
#              random=~vsr(Name,Gu=A),
#              rcov=~units,
#              gTerm = "u:Name",
#              M=marks, data=DT2)
# plot(ans2$scores[,1])
}