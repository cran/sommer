\name{AR1.mat}
\alias{AR1.mat}

\title{
Autocorrelation matrix of order 1.
}
\description{
Creates an autocorrelation matrix of order one with parameters specified. 
}
\usage{
AR1.mat(ar, dimo)
}

\arguments{
  \item{ar}{rho value for the matrix.}
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
R1 <- AR1.mat(.25,4)
image(R1)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
