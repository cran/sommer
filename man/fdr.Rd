\name{fdr}
\alias{fdr}
\title{False Discovery Rate calculation}
\description{
This function calculates the FDR rate of 0.05 and is calculated using the qvalue package (Storey and Tibshirani 2003). The p-value corresponding to a q-value of 0.05 is determined by interpolation. The input is just a vector of p-values.
}
\usage{
fdr(p, fdr.level = 0.05, type = "-log10")
}

\arguments{
  \item{p}{a vector of -log(p-values) or p-values}
  \item{fdr.level}{the level of false discovery rate desired}
  \item{type}{the type of p-values we are entering, either "-log10" or "raw" p-values}
}
\value{
\describe{

\item{$fdr.value}{ a scalar value indicating where the FDR line should be drawn for your raw p-values}
}
}
\references{
Storey and Tibshirani. 2003. Statistical significance for genome-wide studies. PNAS 100:9440-9445
}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
pp <- rnorm(100,15,2) # -log10(p-values) obtained by GWAS
fdr(pp)
}