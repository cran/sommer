\name{us}
\alias{us}
\title{us functionality}
\description{
  This function creates internally the incidence matrices for an unstructure model when using mmer2.
}
\usage{
  us(x)
}

\arguments{
  \item{x}{a column of a dataframe.}
}
\value{
  \describe{
    
    \item{$dd}{ an incidence musrix with columns for the specific levels.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitusive traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{

data(example)
mix <- mmer2(Yield~1, 
              random=~at(Env):Block + us(Env):Name, 
              rcov=~at(Env):units, data=example)
summary(mix)
## compare with asreml, you'll get the same result
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}