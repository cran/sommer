\name{score.calc}
\alias{score.calc}
\title{Score calculation for markers}
\description{
This function is a wrapper from the rrBLUP package to be used when a mixed model including markers to perform GWAS is specified and once the variance components have been estimated the fixed effects are obtained as B= (X'V-X)-X'V-y and the score calculation is obtained with the F statistic as F = Beta^2 / Var(Beta) where Var(Beta) = SSe/(n-p) * [XH-X']-, and quantile value for the beta distribution is calculated as q = (n-p) / (n-p + 1 * F)  which once obtained, the -log10 for such value is the score value.
}
\usage{
score.calc(M, y, ZO, Hinv, min.MAF, X2, P3D=TRUE, method="AI", 
           iters=50, R=NULL, REML=TRUE, draw=TRUE)
}

\arguments{
  \item{M}{marker matrix}
  \item{y}{response variable}
  \item{ZO}{incidence matrix of random effects}
  \item{Hinv}{inverse of the phenotypic variance matrix}
  \item{min.MAF}{minimum minor allele frequency}
  \item{X2}{incidence matrix X as full rank from eigen decomposition}
  \item{P3D}{When P3D=TRUE, variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. The package currently is supported by 3 algorithms; "EMMA" efficient mixed model association (Kang et al. 2008), "AI" average information (Gilmour et al. 1995; Lee et al. 2015), and "EM" expectation maximization (Searle 1993; Bernardo 2010). The default method is average information "AI" because of its ability to handle multiple random effects and its greater speed compared to "EM" which can handle multiple random effects but it is much slower and is not based on REML. On the other hand the "EMMA" method has been implemented but the user needs to keep in mind that such algorithm is based on the theory of ridge regression and the ridge parameter allows the estimation of only a single variance component other than the error.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM algorithm is selected. There is no rule of tumb for the number of iterations. The default value is 50 iterations or EM steps, but could take less or much longer than that.}
  \item{R}{a matrix for variance-covariance structures for the residuals, i.e. for longitudinal data. if not passed is assumed an identity matrix.}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The dafault is TRUE.}

}

\value{
\describe{
\item{$score}{ a vector with the -log10(p-values) for the marker effects in the trait under study}
}
}
\examples{
# it works internally in the \code{\link{mixed.solver}} function
}