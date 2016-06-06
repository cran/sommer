\name{my.colors}
\alias{my.colors}
\docType{data}
\title{
All typical colors in R easy to access.
}
\description{
This dataset is just a vector of the different colors in R from the pdf available at CRAN. Just to make easier the access to them.

}
\usage{data("my.colors")}
\format{
  The format is:
 chr "my.colors"
}
\source{
This data is from CRAN.
}
\references{
Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. https://cran.rstudio.com/web/packages/sommer/ 
}
\examples{
####=========================================####
####=========================================####
data(my.colors)
set.seed(1234)
palette <- sample(my.colors, 16) # sample some
palette2 <- transp(palette) # make them transparent
ma <- matrix(1:16,4,4)
layout(matrix(1:2,1,2))
image(ma, col=palette)
image(ma, col=palette2)
layout(matrix(1,1,1))
####=========================================####
####=========================================####
}
\keyword{datasets}
\seealso{The core function of the package \code{\link{mmer}}}
