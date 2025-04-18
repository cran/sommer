\name{DT_legendre}
\alias{DT_legendre}
\docType{data}
\title{
Simulated data for random regression
}
\description{
A data frame with 4 columns; SUBJECT, X, Xf and Y to show how to use the Legendre polynomials in the mmes function using a numeric variable X and a response variable Y. 

}
\usage{data("DT_legendre")}
\format{
  The format is:
 chr "DT_legendre"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
 
}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
# you need to install the orthopolynom library to do random regression models
# library(orthopolynom)
# data(DT_legendre)
# DT <- DT_legendre
# mRR2<-mmes(Y~ 1 + Xf
#            , random=~ vsm(usm(leg(X,1)),ism(SUBJECT))
#            , rcov=~units
#            , data=DT)
# summary(mRR2)$varcomp

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmes}} }
