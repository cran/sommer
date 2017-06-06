\name{at}
\alias{at}
\title{at functionality}
\description{
  This function creates internally the incidence matrices for specific locations when using mmer2.
}
\usage{
  at(x, levs)
}

\arguments{
  \item{x}{a column of a dataframe.}
  \item{levs}{levels of the columns to keep after using the model.matrix() function.}
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
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
head(hybrid2)
modFD <- mmer2(Yield~1, random=~at(Location,c("3","4")):GCA2, data=hybrid2)
summary(modFD)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }