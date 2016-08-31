\name{ARMA.mat}
\alias{ARMA.mat}

\title{
Autorregressive moving average matrix
}
\description{
Creates an autorregressive moving average matrix with parameters specified. 
}
\usage{
ARMA.mat(ar, lam, dimo)
}

\arguments{
  \item{ar}{rho value for the matrix.}
  \item{lam}{lambda value of the matrix.}
  \item{dimo}{dimensions of the square matrix.}
}
\details{
This function was created to be used in the R structures of the model. Specially useful for spatial variation trends in the field. The values of the parameters will be optimized when passing initial matrices to the \code{\link{mmer}} function.

}
\value{
If everything is defined correctly the function returns:
\describe{
\item{$nn}{the correlation matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}

\examples{
R1 <- ARMA.mat(.25,.2, 4)
image(R1)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
