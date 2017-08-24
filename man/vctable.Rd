\name{vctable}
\alias{vctable}
\title{vctable functionality}
\description{
  This function creates a table of random effects for a multi trait model to be used internally for the mmer and mmer2 function to map the variance components.
}
\usage{
  vctable(ts,nz,nr)
}

\arguments{
  \item{ts}{number of traits.}
  \item{nz}{number of random effects.}
  \item{nr}{number of residual structures.}
}
\value{
  \describe{
    
    \item{all3}{ table of variance-covariance components for each trait.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{
vctable(2,3,1)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }