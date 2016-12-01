\name{MMERM}
\alias{MMERM}
\title{
Multivariate mixed model solver to be called inside mmer
}
\description{
This function estimates parameters for a model of the form \deqn{Y=BX+GZ+ E} where \eqn{Y} is the \eqn{n.obs x r.responses} matrix of response variable, \eqn{X} is a \eqn{n.obs x q.fix.effects} known design matrix of fixed effects,  \eqn{Z} is a \eqn{n.obs x l.ran.effects} known design matrix of random effects, \eqn{B} is \eqn{q.fix.effects x n.reponses} matrix of fixed effects coefficients and  \eqn{G} and \eqn{E} are independent matrix variate variables with \eqn{N_{d x l}(0, V_G, K)} and \eqn{N_{d x n}(0, V_E, I_n)} correspondingly. It also produces the BLUPs for the random effects G and some other statistics.
}
\usage{
MMERM(Y, X=NULL, Z=NULL, W=NULL, method="NR", tolpar = 1e-06, tolparinv = 1e-06,
      draw=TRUE, REML=TRUE, silent=FALSE, iters=15, init=NULL, che=TRUE,
      EIGEND=FALSE, forced=NULL, P3D=TRUE, models="additive", ploidy=2, 
      min.MAF=0.05, gwas.plots=TRUE, map=NULL,manh.col=NULL,fdr.level=0.05,
      constraint=TRUE, IMP=TRUE)
}
\arguments{
   \item{Y}{a numeric vector for the response variable or a data frame of multiple responses (it can run multivariate and also runs univariate in parallel by using the 'n.cores' argument). The multivariate version is only enabled for a single random effect.}
  \item{X}{an incidence matrix for fixed effects related to environmental effects or experimental design. This has to be provided as a matrix, NOT in a list structure.}
  \item{Z}{incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. THIS NEEDS TO BE PROVIDED AS A 2-LEVEL LIST STRUCTURE. For example:
  
  '
  
  ETA <- list(
  
  A=list(Z=Z1, K=K1),
  
  B=list(Z=Z2, K=K2),
  
  C=list(Z=Z3, K=K3)
  
  ) 
  
  
  '
  
  makes a 2 level list for 3 the random effects A, B and C, stored in a variable we call ETA. The general idea is that each random effect is a list, i.e. A=list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix for the random effect, if K is not provided is assumed an identity matrix conferring independence. 
  
  '
  
  PLEASE remember to use the names Z and K FOR ALL RANDOM EFFECTS when you provide your matrices, that's the only way the program distinguishes between a Z or a K matrix. 
  
  '
  
  To provide extra detail, I'll rephrase it; when moving to situations of more than one random effect, you need to build a list for each random effect, and at the end everything gets joined in a list as well (BGLR type of format). Is called a 2-level list, i.e. A=list(Z=Z1, K=K1) and B=list(Z=Z2, K=K2) refers to 2 random effects and they should be put together in a list:
  
  '
  
  ETA <- list( A=list(Z=Z1, K=K1), B=list(Z=Z1, K=K1) )
  
  '
  
  Now you can fit your model as:
  
  '
  
  mod1 <- mmer(y=y, Z=ETA) 
  
  '
  
  You can see the examples at the bottom to have a clearer idea how to fit your models.
  
  }
  
  \item{W}{an incidence matrix for extra fixed effects and only to be used if GWAS is desired and markers will be treated as fixed effects according to Yu et al. (2006) for diploids, and Rosyara et al (2016) for polyploids. Theoretically X and W are both fixed effects, but they are separated to perform GWAS in a model y = Xb + Zu + Wg, allowing the program to recognize the markers from other fixed factors such as environmental factors. This has to be provided as a matrix same than X.
  
  Performs genome-wide association analysis based on the mixed model (Yu et al. 2006):

\deqn{y = X \beta + Z g + W \tau + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The terminology "P3D" (population parameters previously determined) was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the p.adjust function included in the stats package.}

  \item{method}{multivariate method to apply. Currently only "EMMA-M" is available for a single random effect.}
  \item{tolpar}{tolerance parameter for convergence}
  \item{tolparinv}{tolerance parameter for matrix inverse}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM is performed. There is no rule of tumb for the number of iterations. The default value is 100 iterations or EM steps.}
      \item{init}{vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.}
      \item{che}{a TRUE/FALSE value indicating if list structure provided by the user is correct to fix it. The default is TRUE but is turned off to FALSE within the mmer function which would imply a double check.}
      \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this methos are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
      \item{forced}{a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
    \item{P3D}{when the user performs GWAS, P3D=TRUE means that the variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}
    \item{models}{The model to be used in GWAS. The default is the additive model which applies for diploids and polyploids but the model can be a vector with all possible models, i.e. "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref" models are supported for polyploids based on Rosyara (2016).}
      \item{ploidy}{A numeric value indicating the ploidy level of the organism. The default is 2 which means diploid but higher ploidy levels are supported. This should only be modified if you are performing GWAS in polyploids.}
  \item{min.MAF}{when the user performs GWAS min.MAF is a scalar value between 0-1 indicating what is theminor allele frequency to be allowed for a marker during a GWAS analysis when providing the matrix of markers W. In general is known that results for markers with alleles with MAF < 0.05 are not reliable unless sample size is big enough.}
  \item{gwas.plots}{a TRUE/FALSE statement indicating if the GWAS and qq plot should be drawn or not. The default is TRUE but you may want to avoid it.}
  \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' not neccesarily with same dimensions that markers. The program will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
  \item{manh.col}{a vector with colors desired for the manhattan plot. Default are cadetblue and red alternated.}
  \item{fdr.level}{a level of FDR to calculate and plot the line in the GWAS plot. Default is fdr.level=0.05}
  \item{constraint}{a TRUE/FALSE value indicating if the function should apply the boundary constraint indicating that variance components that are zero should be removed from the analysis and variance components recalculated.}
  \item{IMP}{a TRUE/FALSE statement if the function should impute the Y matrix/dataframe with the median value or should get rid of missing values for the estimation of variance components. Only for multivariate mixed models.}
    
}

\value{
\item{Vg }{Estimate of \eqn{V_G}}
\item{Ve}{Estimate of \eqn{V_E}}
\item{Bhat}{BLUEs for \eqn{B}}
\item{Gpred}{BLUPs for \eqn{G}}
\item{XsqtestB}{\eqn{\chi^2} test statistics for testing whether the fixed effect coefficients are equal to zero.}
\item{pvalB}{pvalues obtained from large sample theory for the fixed effects. We report the pvalues adjusted by the "padjust" function for all fixed effect coefficients.}
\item{XsqtestG}{\eqn{\chi^2} test statistic values for testing whether the BLUPs are equal to zero.}
\item{pvalG}{pvalues obtained from large sample theory for the BLUPs. We report the pvalues adjusted by the "padjust" function.}
\item{varGhat}{Large sample variance for BLUPs.}
\item{varBhat}{Large sample variance for the elements of B.}
\item{PEVGhat}{Prediction error variance estimates for the BLUPs.}
}



\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
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
#ans.A <- MMERM(Y=Y, Z=ETA.A)
#ans.A$var.comp
}
