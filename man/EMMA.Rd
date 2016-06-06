\name{EMMA}
\alias{EMMA}
\title{Efficient Mixed Model Association Algorithm}
\description{
This function is used internally in the function \code{\link{mmer}} when ONLY ONE variance component other than the error needs to be estimated through the use of the efficient mixed model association (EMMA) algorithm.
}
\usage{
EMMA(y, X=NULL, Z=NULL, K=NULL, REML=TRUE, silent=FALSE,EIGEND=FALSE)
}
\arguments{
  \item{y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects.}
  \item{Z}{an incidence matrix for the random effect fitted other than the error variance component.}
  \item{K}{a var-cov matrix for the random effect fitted in Z.}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
  \item{EIGEND}{a TRUE/FALSE value indicating if the function should perform the eigen decomposition of the K matrix for square problems to accelerate inversion and estimate parameters faster. Is 2 to 3 times faster than regular EMMA but only works for square problems.}
  
}
\details{
This algorithm is based on Kang et al. (2008), it is based on REML using the ridge regression parameter lambda as the ration of Var(e)/Var(u). This handles models of the form:

.

y = Xb + Zu + e

.

b ~ N[b.hat, 0]  zero variance because is a fixed term

u ~ N[0, K*sigma(u)]  where: K*sigma(u) = G

e ~ N[0, I*sigma(e)]  where: I*sigma(e) = R

y ~ N[Xb, var(Zu+e)]   where;

var(y) = var(Zu+e) = ZGZ+R = V which is the phenotypic variance

.

The function allows the user to specify the incidence matrices for the random effect "u" as Z and its variance-covariance matrix as K, only one Z and K for one random effect.

.

The likelihood function optimized in this algorithm is:

.

logL = (n - p) * log(sum(eta^2/{ lambda + delta})) + sum(log(lambda + delta))

.

where: 
(n-p)  refers to the degrees of freedom
lambda are the eigenvalues mentioned by Kang et al.(2008) 
delta is the REML estimator of the ridge parameter

.

The algorithm can be summarized in the next steps:

.

1) provide initial value for the ridge parameter

2) estimate S = I - X(X'X)-X'

3) obtain the phenotypic variance V = ZKZ' + delta.prov*I 

4) perform an eigen decomposition of SVS 

5) create "lambda"" as the eigenvalues of SVS and "U"" as the eigenvectors

6) estimate eta=U'y 

7) optimize the likelihood shown above providing "eta", "lambdas" and optimize with respect to "delta" which is the ridge parameter and contains Ve/Vu


}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{

\item{$Vu}{ a scalar value for the variance component estimated}
\item{$Ve}{ a scalar value for the error variance estimated}
\item{$V.inv}{ a matrix with the inverse of the phenotypic variance V = ZGZ+R, V^-1}
\item{$u.hat}{ a vector with BLUPs for random effects}
\item{$Var.u.hat}{ a vector with variances for BLUPs}
\item{$PEV.u.hat}{ a vector with predicted error variance for BLUPs}
\item{$beta.hat}{ a vector for BLUEs of fixed effects}
\item{$Var.beta.hat}{ a vector with variances for BLUEs}
\item{$X}{ incidence matrix for fixed effects, if not passed is assumed to only include the intercept}
\item{$Z}{ incidence matrix for random effects, if not passed is assumed to be a diagonal matrix}
\item{$K}{ the var-cov matrix for the random effect fitted in Z}
\item{$ll}{ the log-likelihood value for obtained when optimizing the likelihood function when using ML or REML}

}
}
\references{

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### breeding values with 1 variance component
#### using EMMA algorithm
####=========================================####

####=========================================####
#### random population of 200 lines with 1000 markers
####=========================================####
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- ifelse(runif(1000)<0.5,-1,1)
}
####=========================================####
#### Simulate random phenotypes
####=========================================####
u <- rnorm(1000)
g <- as.vector(crossprod(t(M),u))
h2 <- 0.5  #heritability
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

ans <- EMMA(y=y, K=A.mat(M))
}
\seealso{The core function of the package \code{\link{mmer}}}
