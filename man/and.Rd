\name{and}
\alias{and}
\title{and functionality}
\description{
  overlay (add) r times the design matrix for model term t to the existing design matrix. Specifically, if the model up to this point has p effects and t has a effects, the a columns of the design matrix for t are multiplied by the scalar r (default value 1.0). This can be used to force a correlation of 1 between two terms as in a diallel analysis.
}
\usage{
  and(x)
}

\arguments{
  \item{x}{the random effect to overlay.}
}
\value{
  \describe{
    \item{$x}{ the random effect and the signal to mmer2 to do the overlay.}
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
data(HDdata)
head(HDdata)
HDdata$female <- as.factor(HDdata$female)
HDdata$male <- as.factor(HDdata$male)
HDdata$geno <- as.factor(HDdata$geno)
#### model using overlay
modh <- mmer2(sugar~1, random=~female + and(male) + geno, 
              data=HDdata)
summary(modh)

##################################################################
#### model using overlay [and(.)] and covariance structures [g(.)]
##################################################################
# A <- diag(7); A[1,2] <- 0.5; A[2,1] <- 0.5 # fake covariance structure
# colnames(A) <- as.character(1:7); rownames(A) <- colnames(A);A
# 
# modh2 <- mmer2(sugar~1, random=~g(female) + and(g(male)) + geno, 
#               G=list(female=A, male=A),data=HDdata)
# summary(modh2)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }