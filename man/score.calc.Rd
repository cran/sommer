\name{score.calc}
\alias{score.calc}
\title{Score calculation for markers}
\description{
This function is a wrapper from the rrBLUP package to be used when a mixed model including markers to perform GWAS is specified and once the variance components have been estimated the fixed effects are obtained as B= (X'V-X)-X'V-y and the score calculation is obtained with the F statistic as F = Beta^2 / Var(Beta) where Var(Beta) = SSe/(n-p) * [XH-X']-, and quantile value for the beta distribution is calculated as q = (n-p) / (n-p + 1 * F)  which once obtained, the -log10 for such value is the score value.
}
\usage{
score.calc(marks,y,Z,X,K,ZZ,M,Hinv,ploidy,model,min.MAF,
           max.geno.freq,silent=FALSE,P3D=TRUE)
}

\arguments{
  \item{marks}{marker names}
  \item{y}{response variable}
  \item{Z}{incidence matrix of random effects}
  \item{X}{incidence matrix X as full rank from eigen decomposition}
  \item{K}{covariance structure for random effects}
  \item{ZZ}{incidence matrix of random effects}
  \item{M}{marker matrix}
  \item{Hinv}{inverse of the phenotypic variance matrix}
  \item{ploidy}{numeric value of ploidy level, i.e. 2}
  \item{model}{model for GWAS}
  \item{min.MAF}{minimum minor allele frequency}
  \item{max.geno.freq}{1 - min.MAF}
  \item{silent}{a TRUE/FALSE value indicating if the progress bar should be drawn or not}
  \item{P3D}{when the user performs GWAS, P3D=TRUE means that the variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}

}

\value{
\describe{
\item{$score}{ a vector with the -log10(p-values) for the marker effects in the trait under study}
}
}
\examples{
# it works internally in the \code{\link{mmer}} function
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}