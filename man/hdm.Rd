\name{hdm}
\alias{hdm}
\title{Half Diallel Matrix}
\description{

THIS FUNCTION HAS BEEN OVERPASSED BY THE 'overlay' FUNCTION WHICH IS MORE FLEXIBLE AND FASTER.

This function uses a dataframe with 2 columns named "female" and "male" in numeric format and creates an incidence matrix for a single explanatory variable corresponding to the GCA effect. The resulting incidence matrix can be used in the mmer function as a 'Z' argument in the 2-level list argument for random effects.
}
\usage{
hdm(dat)
}

\arguments{
  \item{dat}{a dataframe with 2 columns named 'female' and 'male' with numeric or factor values indicating the male or female used to produce such hybrid.}
}
\value{
\describe{

\item{$Z}{ an incidence matrix with as many columns as parents in the dataframe indicating with ones the parents used for a particular hybrid (in rows).}
}
}
\references{
Fikret Isik. 2009. Analysis of Diallel Mating Designs. North Carolina State University, Raleigh, USA.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\author{
Giovanny Covarrubias-Pazaran
}

\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}