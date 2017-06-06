\name{phase.F1}
\alias{phase.F1}
\title{
Phasing F1 (CP) data in biparental populations
}
\description{
This function was designed to provide a useful tool for building dense genetic maps in F1 (sometimes referred as CP populations) crosses. This tool alone doesn't build genetic maps but is an important intermediate step researchers will find difficult to deal with.

It takes CP data in a preliminary order (produced by JoinMap4.1 or OneMap) and produces parental maps phased including hk markers translated to parental format without missing the information from hk markers (JoinMap usually looses the information from hk alleles from hk markers, so much missing data makes difficult the bin mapping).

Ideally mapping in F1s would include the following steps:

1) Do initial map with any algorithm able to deal with CP crosses (i.e. JoinMap v4.1 ML algorithm can do maps with all markers lmxll, nnxnp, hkxhk, efxeg, abxcd simultaneously)

\strong{2) Develop parental maps (this is the step this function performs, including phasing and imputation of hk markers)}

3) Do bin mapping in the parental maps (ASMap in R is a good choice for doing the bin maps)

4) Do an integrated map using the parental bin maps (LP merge is a good option)

}
\usage{
phase.F1(genos, silent=FALSE, start=7)
}
\arguments{
  \item{genos}{a dataframe with joinmap format as shown in the \code{\link{F1geno}} data which is explained with more detail in the example below.}
  \item{silent}{a TRUE/FALSE statement indicating if the program should print the progrss bar as the function works.}
  \item{start}{the marker where the program should start to look for the first hk marker.}
  
}
\value{
\describe{

\item{$ff}{ a list with 6 data frames, consensus and parental maps in joinmap and Rqtl format.}

}
}

\author{
Brandon Schlautman
}

\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
data(F1geno)
F1geno[1:10,1:10]
#maps <- phase.F1(F1geno)
#maps$maternal[1:10,1:10]
#heat <- apply(maps$maternal[,-c(1:6)],2,function(x){as.numeric(as.factor(x))})
#heat[1:10,1:10]
#image(heat)
####=========================================####
#### use variance to remove bad markers
####=========================================####
#plot(apply(heat,1,var))
#jjj <- which(apply(heat,1,var) > .5)
#image(heat[-jjj,])
####=========================================####
#### do bin mapping using ASMap
####=========================================####
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}