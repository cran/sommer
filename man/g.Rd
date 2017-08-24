\name{g}
\alias{g}
\title{g functionality}
\description{
  This function allows mmer2 to recognize the use of a variance-covariance structure for a random effect.
}
\usage{
  g(x)
}

\arguments{
  \item{x}{random effect with variance-covariance structure to add}
}
\value{
  \describe{
    
    \item{$x}{ returns the random effect and mmer2 can include the G matrix.}
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
#### call the data
data(CPdata)
#### create the variance-covariance matrix 
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)
mix1 <- mmer2(Yield~1,random=~g(id), G=list(id=A), data=CPpheno)
summary(mix1)

}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}