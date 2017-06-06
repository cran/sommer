\name{eigenGWAS}
\alias{eigenGWAS}
\title{Unraveling selection signatures with eigenGWAS}
\description{
This function performs eigenGWAS based on the paper proposed by Chen et al. (2016) to avoid expensive computations with high number of markers. For mathematical details of the derivation please refer to the original publication below in the references.

eigenGWAS is designed to find ancestry informative markers (AIM) and loci under selection. The phenotype of eigenGWAS is an eigenvector generated from marker data.

}
\usage{
eigenGWAS(markers, eivec=1, map=NULL)
}

\arguments{
  \item{markers}{a matrix of numeric markers}
  \item{eivec}{eigen vector to be used as a response}
  \item{map}{optional argument to provide a map to sort the signature. Obligatory column names of the map are "Chrom" and "Position".}
}
\value{
\describe{

\item{$result}{ list with signature and Fst metric.}

}
}
\references{
Chen, G.B., S.H. Lee, ZX Zhu, B Benyamin, MM Robinson, EigenGWAS: finding loci under selection through genome-wide association studies of eigenvectors in structured populations. Heredity 2016.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
data(CPdata)
markers <- CPdata$geno
markers[1:5,1:5]
#res <- eigenGWAS(markers)
### ===================== ###
### use map information
### ===================== ###
#my.map <- CPdata$map
#head(my.map)
#res <- eigenGWAS(markers, map = my.map)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}