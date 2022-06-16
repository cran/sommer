\name{csc}
\alias{csc}
\title{customized covariance structure}
\description{
  \code{csc} creates a customized covariance structure for specific levels of the random effect to be used with the \code{\link{mmec}} solver.
}
\usage{
  csc(x, mm, thetaC, theta)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
  \item{mm}{customized variance-covariance structure for the levels of the random effect.}
  \item{thetaC}{an optional matrix for constraints in the variance components.}
  \item{theta}{an optional matrix for initial values of the variance components.}
}
\value{
  \describe{
    
    \item{$res}{ a list with the provided vector and the variance covariance structure expected for the levels of the random effect.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
x <- as.factor(c(1:5,1:5,1:5));x
csc(x,matrix(1,5,5))
}
\seealso{The function \code{\link{vsc}} to know how to use \code{csc} in the \code{\link{mmec}} solver.}