\name{vctable.help}
\alias{vctable.help}
\title{vctable.help functionality}
\description{
  This function extracts the structure for multi trait models for the random and rcov formulas of mmer2.
}
\usage{
  vctable.help(random=NULL,rcov=NULL, data)
}

\arguments{
  \item{random}{formula for random terms.}
  \item{rcov}{formula for residual terms.}
  \item{data}{dataset.}
}
\value{
  \describe{
    
    \item{res}{trait structures for each random and residual term.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }