\name{DT_technow}
\alias{DT_technow}
\alias{Ad_technow}
\alias{Af_technow}
\alias{Md_technow}
\alias{Mf_technow}
\docType{data}
\title{
Genotypic and Phenotypic data from single cross hybrids (Technow et al.,2014)
}
\description{
This dataset contains phenotpic data for 2 traits measured in 1254 single cross hybrids coming from the cross of Flint x Dent heterotic groups. In addition contains the genotipic data (35,478 markers) for each of the 123 Dent lines and 86 Flint lines. The purpose of this data is to demosntrate the prediction of unrealized crosses (9324 unrealized crosses, 1254 evaluated, total 10578 single crosses). We have added the additive relationship matrix (A) but can be easily obtained using the A.mat function on the marker data. Please if using this data for your own research cite Technow et al. (2014) publication (see References).

}
\usage{data("DT_technow")}
\format{
  The format is:
 chr "DT_technow"
}
\source{
This data was extracted from Technow et al. (2014).
}
\references{
If using this data for your own research please cite:

Technow et al. 2014. Genome properties and prospects of genomic predictions of hybrid performance in a Breeding program of maize. Genetics 197:1343-1355.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(DT_technow)
DT <- DT_technow
Md <- Md_technow
Mf <- Mf_technow
# Md <- (Md*2) - 1
# Mf <- (Mf*2) - 1
# Ad <- A.mat(Md)
# Af <- A.mat(Mf)
# ###=========================================####
# ###=========================================####
# ans2 <- mmes(GY~1,
#              random=~vsm(ism(dent),Gu=Ad) + vsm(ism(flint),Gu=Af),
#              rcov=~units,
#              data=DT)
# summary(ans2)$varcomp
# 
# Adi <- solve(Ad + diag(1e-4,ncol(Ad),ncol(Ad)))
# Adi <- as(as(as( Adi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Adi, 'inverse')=TRUE 
# Afi <- solve(Af + diag(1e-4,ncol(Af),ncol(Af)))
# Afi <- as(as(as( Afi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Afi, 'inverse')=TRUE 
# ####====================####
# #### multivariate model ####
# ####     2 traits       ####
# ####====================####
# head(DT)
# DT2 <- stackTrait(DT, traits = c("GY","GM"))
# head(DT2$long)
# 
# M <- rbind(Md,Mf)
# A <- A.mat(M)
# Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, 'inverse')=TRUE 
# DT2$long=DT2$long[with(DT2$long, order(trait)), ]
# 
# ans3 <- mmes(valueS~trait, henderson=TRUE,
#              random=~vsm(usm(trait),ism(overlay(dent,flint)),Gu=Ai),
#              rcov=~ vsm(dsm(trait), ism(units)),
#              data=DT2$long)
# summary(ans3)
# cov2cor(ans3$theta[[1]])
}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmes}} }
