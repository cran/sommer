\name{rrm}
\alias{rrm}
\title{reduced rank covariance structure}
\description{
  \code{rrm} creates a reduced rank factor analytic covariance structure by selecting the n vectors of the L matrix of the Cholesky decomposition or the U vectors of the SVD decomposition (loadings or latent covariates) to create a new incidence matrix of latent covariates that can be used with the \code{\link{mmes}} solver to fit random regressions on the latent covariates.
}
\usage{
  rrm(x=NULL, H=NULL, nPC=2, returnGamma=FALSE, cholD=TRUE)
}

\arguments{
  \item{x}{vector of the dataset containing the variable to be used to form the incidence matrix.}
  \item{H}{two-way table of identifiers (rows; e.g., genotypes) by features (columns; e.g., environments) effects. Row names and column names are required. No missing data is allowed.}
  \item{nPC}{number of principal components to keep from the loadings matrix.}
  \item{returnGamma}{a TRUE/FALSE argument specifying if the function should return the matrix of loadings used to build the incidence matrix for the model. The default is FALSE so it returns only the incidence matrix.}
  \item{cholD}{a TRUE/FALSE argument specifying if the Cholesky decomposition should be calculated or the singular value decomposition should be used instead.}
  
}
\details{
This implementation of a version of the reduced rank factor analytic models uses the so-called principal component (PC) models (Meyer, 2009) which assumes specific effects (psi) are equal to 0. The model is as follows:

y = Xb + Zu + e

where the variance of u ~ MVN(0, Sigma)

Sigma = (Gamma_t Gamma) + Psi

\strong{Extended factor analytic model}:

y = Xb + Z(I Gamma)c + Zs + e = Xb + Z*c + Zs + e

where \code{y} is the response variable, \code{X} and \code{Z} are incidence matrices for fixed and random effects respectively, \code{I} is a diagonal matrix, \code{Gamma} are the factor loadings for \code{c} common factor scores, and \code{s} are the specific effects, \code{e} is the vector of residuals.   

\strong{Reduced rank model:}

y = Xb + Z(I Gamma)c + e = Xb + Z*c + e

which is equal to the one above but assumes specific effects = 0. 

\strong{The algorithm in rrm is the following:}

1) uses a wide-format table of timevar (m columns) by idvar (q rows) named H to form the initial variance-covariance matrix (Sigma) which is calculated as Sigma = H'H of dimensions m x m (column dimensions, e.g., environments x environments). 

2) The Sigma matrix is then center and scaled. 

3) A Cholesky (L matrix) or SVD decomposition (U D V') is performed in the Sigma matrix.

4) n vectors from L (when Cholesky is used) or U sqrt(D) (when SVD is used) are kept to form Gamma. Gamma = L[,1:nPc] or Gamma = U[,1:nPC]. These are the so-called loadings (L for all loadings, Gamma for the subset of loadings).

4) Gamma is used to form a new incidence matrix as Z* = Z Gamma

5) This matrix is later used for the REML machinery to be used with the usc (unstructured) or dsm (diagonal) structures to estimate variance components and factor scores. The resulting BLUPs from the mixed model are the optimized factor scores. Pretty much as a random regression over latent covariates.

This implementation does not update the loadings (latent covariates) during the REML process, only estimates the REML factor scores for fixed loadings. This is different to other software (e.g., asreml) where the loadings are updated during the REML process as well. 

BLUPs for genotypes in all locations can be recovered as:

u = Gamma * u_scores

The resulting loadings (Gamma) and factor scores can be thought as an equivalent to the classical factor analysis.

}
\value{
  \describe{
    
    \item{$Z}{ a incidence matrix Z* = Z Gamma which is the original incidence matrix for the timevar multiplied by the loadings.}
    
    \item{$Gamma}{ a matrix of loadings or latent covariates.}
        
    \item{$Sigma}{ the covariance matrix used to calculate Gamma.}
    
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744
  
  Meyer K (2009) Factor analytic models for genotype by environment type problems and structured covariance matrices. Genetics Selection Evolution, 41:21
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_h2)
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]
head(DT)
indNames <- na.omit(unique(DT$Name))
A <- diag(length(indNames))
rownames(A) <- colnames(A) <- indNames

\donttest{
  
# fit diagonal model first to produce H matrix
ansDG <- mmes(y~Env, henderson=TRUE,
              random=~ vsm(dsm(Env), ism(Name)),
              rcov=~units, nIters = 100,
              data=DT)

H0 <- ansDG$uList$`vsm(dsm(Env), ism(Name))` # GxE table

# reduced rank model
ansFA <- mmes(y~Env, henderson=TRUE,
              random=~vsm( usm(rrm(Env, H = H0, nPC = 3)) , ism(Name)) + # rr
                vsm(dsm(Env), ism(Name)), # diag
              rcov=~units,
              # we recommend giving more iterations to these models
              nIters = 100,
              data=DT)

vcFA <- ansFA$theta[[1]]
vcDG <- ansFA$theta[[2]]

loadings=with(DT, rrm(Env, nPC = 3, H = H0, returnGamma = TRUE) )$Gamma
scores <- ansFA$uList[[1]]

vcUS <- loadings \%*\% vcFA \%*\% t(loadings)
G <- vcUS + vcDG
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)

uFA <- scores \%*\% t(loadings)
uDG <- ansFA$uList[[2]]
u <- uFA + uDG
  
}


}
\seealso{The function \code{\link{vsm}} to know how to use \code{rrm} in the \code{\link{mmes}} solver.}