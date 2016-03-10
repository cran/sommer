\name{AI}
\alias{AI}
\title{Average Information Algorithm}
\description{
This function is used internally in the function \code{\link{mmer}} when MORE than 1 variance component needs to be estimated through the use of the average information (AI) algorithm.
}
\usage{
AI(y, X=NULL, ZETA=NULL, R=NULL, draw=TRUE, REML=TRUE, silent=FALSE, 
   iters=50, constraint=TRUE, init=NULL, sherman=FALSE, che=TRUE, 
   MTG2=FALSE, Fishers=FALSE)
}
\arguments{
  \item{y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects.}
  \item{ZETA}{an incidence matrix for random effects. This can be for one or more random effects. This NEEDS TO BE PROVIDED AS A LIST STRUCTURE. For example Z=list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=K3)) makes a 2 level list for 3 random effects. The general idea is that each random effect with or without its variance-covariance structure is a list, i.e. list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix. When moving to more than one random effect we need to make several lists that need to be inside another list. What we call a 2-level list, i.e. list(Z=Z1, K=K1) and list(Z=Z2, K=K2) would need to be put in the form; list(list(Z=Z1, K=K1),list(Z=Z1, K=K1)), which as can be seen, is a list of lists (2-level list). }
  \item{R}{a matrix for variance-covariance structures for the residuals, i.e. for longitudinal data. if not passed is assumed an identity matrix.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM is performed. There is no rule of tumb for the number of iterations. The default value is 100 iterations or EM steps.}
    \item{constraint}{a TRUE/FALSE value indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
      \item{init}{vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.}
      \item{sherman}{a TRUE/FALSE value indicating if Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467) should be used when estimating variance components in order to perform faster when a mixed model with no covariance structure using the average information algorithm is fitted. The default is FALSE since this software was designed for unreplicated data (altough can fit models with replicated data but slower than lme4).}
      \item{che}{a TRUE/FALSE value indicating if list structure provided by the user is correct to fix it. The default is TRUE but is turned off to FALSE within the mmer function which would imply a double check.}
      \item{MTG2}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this methos are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
      \item{Fishers}{a TRUE/FALSE value indicating if the program should calculate at the final step and return the inverse of the Fishers Information Matrix.}

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

6) fill the average information matrix (AI) with equation provided in Gilmour et al. (1995)

7) obtain AI.inv by inverting AI (the average information matrix)

8) calculate scores by first derivatives refer as "B" in Gilmour et al. (1995)

9) update the values of variance components by :
k(i+1) = k(i) + [ B(i) * AI.inv ]

10) steps are repeated in a while loop until convergence is reached, the likelihood doesn't increase anymore.


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

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

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
#ans <- AI(y=y, ZETA=ETA)
#ans$var.comp
}
