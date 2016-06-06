\name{design.score}
\alias{design.score}
\title{design score for the model to be tested}
\description{
This function is a wrapper from the GWASpoly package that creates the incidence matrix to be used when calculating the -log10 scores for the markers once the variance components and marker effects are estimated.
}
\usage{
design.score(Mi,model,ploidy,min.MAF,max.geno.freq)
}

\arguments{
  \item{Mi}{marker to be used}
  \item{model}{model to be used for designing the matrix. The default is "additive" which can be extend to any ploidy level but alternative models for polyploids such as "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref" are supported based on Rosyara (2016)}
  \item{ploidy}{A numeric value indicating the ploidy level of the organism. The default is 2 which means diploid but higher ploidy levels are supported.}
  \item{min.MAF}{minimum minor allele frequency}
  \item{max.geno.freq}{1 - min.MAF}
}

\value{
\describe{
\item{$S}{ an incidence matrix for the marker based on the model provided.}
}
}
\examples{
# it works internally in the \code{\link{mmer}} function
}
\seealso{The core function of the package \code{\link{mmer}}}