\name{fdr}
\alias{fdr}
\title{False Discovery Rate calculation}
\description{
This function calculates the false discovery rate (FDR) at a certain value specified by the user. Is a wrapper of the p.adjust function from the stats package. It has been adjusted to facilitate the user to find the FDR value for a vector of values in scale minus log10. It is used internally of the \code{\link{mmer}} function when the argument W is used to perform GWAS.
}
\usage{
fdr(p, fdr.level = 0.05)
}

\arguments{
  \item{p}{a vector of minus log(p values)}
  \item{fdr.level}{the level of false discovery rate desired}
}
\value{
\describe{

\item{$p.ad}{ a vector of new minus log10(p values) adjusted for false discovery rate at a given level specified by the user.}
\item{$fdr.value}{ a scalar value indicating where the FDR line should be drawn for your new adjusted minus log10(p values)}
\item{$p.or}{ a vector with the initial minus log10(p values) provided by the user.}
\item{$fdr.or}{ a scalar value indicating where the FDR line should be drawn for your original minus log10(p values), if the user prefers to plot the original values and wants to draw the FDR line in the original scale.}

}
}
\references{
Benjamini, Y., and Yekutieli, D. 2001. The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165-1188.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
####=========================================####
#### generate your mickey mouse -log10(p-values)
####=========================================####
set.seed(1253)
pp <- abs(rnorm(500,0,3));pp[23:34] <- abs(rnorm(12,0,20))

####=========================================####
#### see how they look like
####=========================================####
plot(pp, col=transp("cadetblue"), pch=20, cex=1.5)

####=========================================####
#### adjust the values for FDR and see how they look like
####=========================================####
new.pp <- fdr(pp)
plot(new.pp$p.ad, col=transp("cadetblue"), pch=20, cex=1.5)

####=========================================####
#### or you may want to plot your original values 
#### with the FDR line instead
####=========================================####
plot(pp, col=transp("cadetblue"), pch=20, cex=1.5)
abline(h=new.pp$fdr.10, lty=3, col="red", lwd=2)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}