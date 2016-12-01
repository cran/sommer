\name{MAI}
\alias{MAI}
\title{Multivariate Average Information Algorithm}
\description{
This function is used internally in the function \code{\link{mmer}} when MORE than 1 variance component needs to be estimated through the use of the average information (MAI) algorithm and the MVM argument is set to TRUE.
}
\usage{
MAI(Y, X=NULL, ZETA=NULL, draw=TRUE, REML=TRUE, silent=FALSE, iters=20, 
    init=NULL, tol=1e-3, che=TRUE, EIGEND=FALSE, forced=NULL, IMP=FALSE)

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
  
  mod1 <- mmer(y=y, Z=ETA) 
  
  '
  
  You can see the examples at the bottom to have a clearer idea how to fit your models. }
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM is performed. There is no rule of tumb for the number of iterations. The default value is 100 iterations or EM steps.}
      \item{init}{vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.}
      \item{tol}{Convergence criteria.  If the change in residual log
    likelihood for one cycle is less than \code{10 x tol} the algorithm
    finishes.  If each component of the change proposed by the
  Newton-Raphson is lower in magnitude than \code{tol} the algorithm
  finishes. Default value is \code{1e-4}.}
      \item{che}{a TRUE/FALSE value indicating if list structure provided by the user is correct to fix it. The default is TRUE but is turned off to FALSE within the mmer function which would imply a double check.}
      \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this methos are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
      \item{forced}{a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
      \item{IMP}{a TRUE/FALSE statement if the function should impute the Y matrix/dataframe with the median value or should get rid of missing values for the estimation of variance components.}

}
\details{
This algorithm is based on Gilmour et al. (1995), it is based on REML. This handles models of the form:

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

The algorithm can be summarized in the next steps:

.

1) provide initial values for the variance components

2) estimate the phenotypic variance matrix V = ZGZ + R

3) obtain Vinv by inverting V 

4) obtain the projection matrix P = Vinv - [Vinv X (X'V-X)- X Vinv]

5) evaluate the logLikelihood as shown above

6) fill the average information matrix (MAI) with equation provided in Gilmour et al. (1995)

7) obtain MAI.inv by inverting MAI (the average information matrix)

8) calculate scores by first derivatives refer as "B" in Gilmour et al. (1995)

9) update the values of variance components by :
k(i+1) = k(i) + [ B(i) * MAI.inv ]

10) steps are repeated in a while loop until convergence is reached, the likelihood doesn't increase anymore.


}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{

\item{$Vu}{a scalar value for the variance component estimated}
\item{$Ve}{a scalar value for the error variance estimated}
\item{$V.inv}{a matrix with the inverse of the phenotypic variance V = ZGZ+R, V^-1}
\item{$u.hat}{a vector with BLUPs for random effects}
\item{$Var.u.hat}{a vector with variances for BLUPs}
\item{$PEV.u.hat}{a vector with predicted error variance for BLUPs}
\item{$beta.hat}{a vector for BLUEs of fixed effects}
\item{$Var.beta.hat}{a vector with variances for BLUEs}
\item{$X}{incidence matrix for fixed effects, if not passed is assumed to only include the intercept}
\item{$ZETA}{ random effect incidence and variace-covariance matrices}
\item{$ll}{ the log-likelihood value for obtained when optimizing the likelihood function when using ML or REML}

}
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Lee et al. 2015. EIGEND: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

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
#ans.A <- MAI(Y=Y, ZETA=ETA.A)
#ans.A$var.comp
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
