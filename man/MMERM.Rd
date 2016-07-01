\name{MMERM}
\alias{MMERM}
\title{
Multivariate mixed model solver
}
\description{
This function estimates parameters for a model of the form \deqn{Y=BX+GZ+ E} where \eqn{Y} is the \eqn{n.obs x r.responses} matrix of response variable, \eqn{X} is a \eqn{n.obs x q.fix.effects} known design matrix of fixed effects,  \eqn{Z} is a \eqn{n.obs x l.ran.effects} known design matrix of random effects, \eqn{B} is \eqn{q.fix.effects x n.reponses} matrix of fixed effects coefficients and  \eqn{G} and \eqn{E} are independent matrix variate variables with \eqn{N_{d x l}(0, V_G, K)} and \eqn{N_{d x n}(0, V_E, I_n)} correspondingly. It also produces the BLUPs for the random effects G and some other statistics.
}
\usage{
MMERM(Y, X=NULL, Z=NULL, method="EMMA", tolpar = 1e-06, 
     tolparinv = 1e-06, che=TRUE, silent=TRUE)
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
  \item{method}{multivariate method to apply. Currently only "EMMA-M" is available for a single random effect.}
  \item{tolpar}{tolerance parameter for convergence}
  \item{tolparinv}{tolerance parameter for matrix inverse}
  \item{che}{a TRUE/FALSE value indicating if list structure provided by the user is correct to fix it. The default is TRUE but is turned off to FALSE within the mmer function which would imply a double check.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
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
