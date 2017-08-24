\name{PolyData}
\alias{PolyData}
\docType{data}
\title{
Genotypic and Phenotypic data for a potato polyploid population
}
\description{
This dataset contains phenotpic data for 18 traits measured in 187 individuals from a potato diversity panel. In addition contains genotypic data for 221 individuals genotyped with 3522 SNP markers. Please if using this data for your own research make sure you cite Rosyara's (2015) publication (see References).

}
\usage{data("PolyData")}
\format{
  The format is:
 chr "PolyData"
}
\source{
This data was extracted from Rosyara (2016).
}
\references{
If using this data for your own research please cite:

Rosyara Umesh R., Walter S. De Jong, David S. Douches, Jeffrey B. Endelman. Software for genome-wide association studies in autopolyploids and its application to potato. The Plant Genome 2015.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(PolyData)
genotypes <- PolyData$PGeno
phenotypes <- PolyData$PPheno

# ####=========================================####
# ####### convert markers to numeric format
# ####=========================================####
# numo <- atcg1234(data=genotypes, ploidy=4); numo[1:5,1:5]; dim(numo)
# 
# ###=========================================####
# ###### plants with both genotypes and phenotypes
# ###=========================================####
# common <- intersect(phenotypes$Name,rownames(numo))
# 
# ###=========================================####
# ### get the markers and phenotypes for such inds
# ###=========================================####
# marks <- numo[common,]; marks[1:5,1:5]
# phenotypes2 <- phenotypes[match(common,phenotypes$Name),];
# phenotypes2 <- as.data.frame(phenotypes2)
# phenotypes2[1:5,]
# 
# ###=========================================####
# ###### Additive relationship matrix, specify ploidy
# ###=========================================####
# A <- A.mat(marks, ploidy=4)
# D <- D.mat(marks, ploidy=4)
# E <- E.mat(marks, ploidy=4)
# ###=========================================####
# ### run as mixed model
# ###=========================================####
# ans <- mmer2(tuber_shape~1, random=~g(Name),
#              G=list(Name=A), data=phenotypes2)
# summary(ans)
# 
# ###=========================================####
# ### run it as GWAS model
# ###=========================================####
# my.map <- PolyData$map
# models <- c("additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
# ans2 <- GWAS2(tuber_shape~1, random=~g(Name), models = "additive",
#              G=list(Name=A), W=marks, data=phenotypes2)
# summary(ans2)
}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
