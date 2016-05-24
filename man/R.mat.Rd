\name{R.mat}
\alias{R.mat}

\title{
Residual relationship matrix
}
\description{
Creates covariance matrices. Creates 3 types of matrices. 1) spatial correlation with no particular direction other than distance, 2) spatial correlation with one direction, 3) spatial correlation with 2 directions. 
}
\usage{
R.mat(dims, rhor=.95, rhoc=.95, method="SIM")
}

\arguments{
  \item{dims}{dimensions of the square matrix.}
\item{rhor}{degree of autocorrelation in the row direction}
\item{rhoc}{degree of autocorrelation in the column direction}
\item{method}{one of the three methods available for creating correlation matrices for random effects related to spetial autocorrelation:
1) `SIM` spatial correlation with no particular direction other than distance, 

2) `AR1` spatial correlation with one direction, 

3) `AR2` spatial correlation with 2 directions. 
}

}
\details{
This function was created to stablish connections between levels of random effects such as row and columns in an experimental design. Hopefully users will find this helpful for adjusting for spatial correlation and avoid the use of correlation structures in the residual part of the model.

}
\value{
If everything is defined correctly the function returns:
\describe{
\item{$nn}{the correlation matrix}
}
}

\references{

Covarrubias-Pazaran G. 2016. Genomic prediction using the R package sommer. Genetics X:yyy-yyy.

}

\examples{
R1 <- R.mat(10, method="SIM")
R2 <- R.mat(10, method="AR2")
layout(matrix(1:2,1,2))
image(R1);image(R2)
}
