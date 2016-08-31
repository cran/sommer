\name{NRR}
\alias{NRR}
\title{Newton-Raphson Algorithm including Residual structures}
\description{
This function is used internally in the function \code{\link{mmer}} when MORE than 1 variance component needs to be estimated through the use of the Newton-Raphson (NR) algorithm allowing the use of residual structures.
}
\usage{
NRR(y,X=NULL,Z=NULL,R=NULL,tolpar=1e-6,tolparinv=1e-6,maxcyc=10,
   draw=TRUE, constraint=TRUE)
}
\arguments{
  \item{y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects.}
  \item{Z}{an incidence matrix for random effects. This can be for one or more random effects. This NEEDS TO BE PROVIDED AS A LIST STRUCTURE. For example Z=list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=K3)) makes a 2 level list for 3 random effects. The general idea is that each random effect with or without its variance-covariance structure is a list, i.e. list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix. When moving to more than one random effect we need to make several lists that need to be inside another list. What we call a 2-level list, i.e. list(Z=Z1, K=K1) and list(Z=Z2, K=K2) would need to be put in the form; list(list(Z=Z1, K=K1),list(Z=Z1, K=K1)), which as can be seen, is a list of lists (2-level list). }
  \item{R}{a two level list including all the R matrices to be included in the analysis. Each element of the two level list is a list with all R matrices to be used. For example, to model spatial variation in a plot with rows(4) and columns(18) you need to create a list with both matrices:
  
  # autocorrelation matrix for the 4 rows, initial value gamma=0.25
  
  \strong{R1 <- AR1.mat(.25,4) }
  
 # autocorrelation matrix for the 18 cols, initial value gamma=0.25
 
  \strong{R2 <- AR1.mat(.25,18)}
  
  then we do the 2-level list as:
  
  \strong{RETA <- list(spatial=list(R1,R2,type=c("AR1","AR1")))}
  
  which can be introduced in the R argument of the function. The idea is that a kronecker product will be taken with R1 and R1, you have to introduce the R matrices that will yield the right dimensions, and you can add as many R matrices as you want. Here we only show using one for spatial effects.
  }
  \item{maxcyc}{a scalar value indicating how many iterations have to be performed if the EM is performed. There is no rule of tumb for the number of iterations. The default value is 100 iterations or EM steps.}
    \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
    \item{constraint}{a TRUE/FALSE value indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
   \item{tolpar}{Convergence criteria.  If the change in residual log
    likelihood for one cycle is less than \code{10 x tol} the algorithm
    finishes.  If each component of the change proposed by the
  Newton-Raphson is lower in magnitude than \code{tol} the algorithm
  finishes. Default value is \code{1e-4}.}
   \item{tolparinv}{Value to be used when the V matrix cannot be inverted so this value will be used to the diagonal of the V matrix to allow inversion.}

}
\details{
This algorithm is based on Tunnicliffe (1989), it is based on REML. This handles models of the form:

.

y = Xb + Zu + e

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

where G is:

.

|K*sigma(gca1).....................0..........................0.........|
    
|.............0.............S*sigma(gca2).....................0.........| = G

|.............0....................0......................W*sigma(sca)..|

.

The likelihood function optimized in this algorithm is:

.

logL = -0.5 * (log( | V | ) + log( | X'VX | ) + y'Py

.

where: | |  refers to the derminant of a matrix

.


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

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

####=========================================####
#### breeding values with 3 variance components
####=========================================####

####=========================================####
## Import phenotypic data on inbred performance
## Full data
####=========================================####
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K # extract the var-cov K

y <- hybrid2$Yield
X1 <- model.matrix(~ Location, data = hybrid2);dim(X1)
Z1 <- model.matrix(~ GCA1 -1, data = hybrid2);dim(Z1)
Z2 <- model.matrix(~ GCA2 -1, data = hybrid2);dim(Z2)
Z3 <- model.matrix(~ SCA -1, data = hybrid2);dim(Z3)

####=========================================####
#### Realized IBS relationships for set of parents 1
####=========================================####
K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)     
####=========================================####
#### Realized IBS relationships for set of parents 2
####=========================================####
K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)     
####=========================================####
#### Realized IBS relationships for cross 
#### (as the Kronecker product of K1 and K2)
####=========================================####
S <- kronecker(K1, K2) ; dim(S)   
rownames(S) <- colnames(S) <- levels(hybrid2$SCA)

ETA <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=S))
####=========================================####
#### run the next line, it was ommited for CRAN time limitations
####=========================================####
#ans <- NRR(y=y, ZETA=ETA)
#ans$var.comp
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
