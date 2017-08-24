\name{eig}
\alias{eig}
\title{eigen decomposition functionality}
\description{
  This function allows mmer2 to detect when the user wants to perform an eigen decomposition of the additive relationship matrix in certain random effect.
}
\usage{
  eig(x)
}

\arguments{
  \item{x}{a column of a dataframe.}
}
\value{
  \describe{
    
    \item{$dd}{ an incidence matrix with columns for the specific levels specified.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{

data(CPdata)
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix 
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)

mix1 <- mmer2(color~1,
              random=~ eig(g(id)), 
              rcov=~units,
              G=list(id=A), 
              data=CPpheno)
summary(mix1)

# mix2 <- mmer2(cbind(Yield,color)~1,
#               random=~ us(trait):eig(g(id)), 
#               rcov=~us(trait):units,
#               G=list(id=A), 
#               data=CPpheno)
# summary(mix2)

}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }