\name{MNR}
\alias{MNR}
\title{Multivariate Newton-Raphson Algorithm}
\description{
This function is used internally in the function \code{\link{mmer}} when MORE than 1 variance component needs to be estimated through the use of the Newton-Raphson (MNR) algorithm for multivariate models (multiple responses).
}
\usage{
MNR(Y,X=NULL,ZETA=NULL,R=NULL,init=NULL,iters=20,tolpar=1e-3,
    tolparinv=1e-6,draw=FALSE,silent=FALSE, constraint=TRUE, 
    EIGEND=FALSE, forced=NULL, IMP=FALSE, complete=TRUE, 
    che=FALSE, restrained=NULL)
}
\arguments{
  \item{Y}{a matrix or data frame of response variables}
  \item{X}{an incidence matrix for fixed effects.}
  \item{ZETA}{incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. \strong{This needs to be provided as a 2-level list structure}. For example:
  
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
 \item{R}{list of R matrices to correct for residual. Internally the program will do the kronecker product of such matrices to create R.}
  \item{init}{vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.}
  \item{iters}{Maximum number of iterations allowed.  Default value is
    15.}
  \item{tolpar}{Convergence criteria.}
  \item{tolparinv}{tolerance parameter for matrix inverse}
    \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
    \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar while working or should not be displayed. The default is FALSE, which means is not silent and will display the progress bar.}
    \item{constraint}{a TRUE/FALSE value indicating if the function should apply the boundary constraint indicating that variance components that are zero should be removed from the analysis and variance components recalculated.}
    \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this method are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
    \item{forced}{a list of values for variance-covariance components to be used if the user wants to force those values.}
    \item{IMP}{a TRUE/FALSE statement if the function should impute the Y matrix/dataframe with the median value or should get rid of missing values for the estimation of variance components.}
   \item{complete}{a TRUE/FALSE statement to indicate if the function should impute the cases where at least for one trait there's an observation.}
   \item{che}{a TRUE/FALSE statement to indicate if the function should check the input parameters from the user.}
   \item{restrained}{a numeric argument specifying which variance-covariance parameters should be restrained.}
}
\details{

Please refer to the \code{\link{sommer}} help page.


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
\item{fish.inv}{inverse of the Fisher's information or average information matrices to obtain variance-covariance of the variance components.}
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
\item{fish.inv}{ If was set to TRUE the Fishers information matrix will be in this slot.}
\item{method}{The method for extimation of variance components specified by the user.}
\item{forced}{a TRUE/FALSE statement indicating if user used the forced argument.}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{restrained}{table of restrained parameters.}

}
\references{

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
data(CPdata)
### look at the data
head(CPpheno)
CPgeno[1:5,1:5]
## fit a model including additive and dominance effects
Y <- CPpheno[,c("color","Yield")]
Za <- diag(dim(Y)[1])
A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- MNR(Y=Y, ZETA=ETA.A)
#ans.A$var.comp
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
