\name{PolyData}
\alias{PolyData}
\docType{data}
\title{
Genotypic and Phenotypic data for a potato polyploid population
}
\description{
This dataset contains phenotpic data for 18 traits measured in 187 individuals from a potato diversity panel. In addition contains genotypic data for 221 individuals genotyped with 3522 SNP markers. Please if using this data for your own research cite Rosyara's (2015) publication (see References).

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

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
data(PolyData)
genotypes <- PolyData$PGeno
phenotypes <- PolyData$PPheno

####=========================================####
####### convert markers to numeric format
####=========================================####
#numo <- atcg1234(data=genotypes, ploidy=4); numo[1:5,1:5]; dim(numo)

####=========================================####
####### plants with both genotypes and phenotypes
####=========================================####
#common <- intersect(phenotypes$Name,rownames(numo))

####=========================================####
#### get the markers and phenotypes for such inds
####=========================================####
#marks <- numo[common,]; marks[1:5,1:5]
#phenotypes2 <- phenotypes[match(common,phenotypes$Name),];
#phenotypes2[1:5,1:5]

####=========================================####
#### response variable
####=========================================####
#yy <- phenotypes2$tuber_shape
#set.seed(1234)
#ww <- sample(1:187,38)
#yy.trn <- yy; yy.trn[ww] <- NA

####=========================================####
####### Additive relationship matrix, specify ploidy
####=========================================####
#K1 <- A.mat(marks, ploidy=4); dim(K1);K1[1:5,1:5]

####=========================================####
#### Incidence matrix for genotypes
####=========================================####
#Z1 <- diag(length(yy))

####=========================================####
#### double check all dimensions
####=========================================####
#dim(Z1); dim(K1); length(yy)
#ETA <- list( list(Z=Z1, K=K1)) # random effects for genotypes

####=========================================####
#### run the genomic selection model
####=========================================####
#ans <- mmer(y=yy.trn, Z=ETA, method="EMMA")
#cor(yy[ww],ans$fitted.y[ww])
#summary(ans)

####=========================================####
#### run it as GWAS model
####=========================================####
#my.map <- PolyData$map
#models <- c("additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
#ans2 <- mmer(y=yy, Z=ETA, W=marks, method="EMMA", 
#              ploidy=4, models=models[1], map=my.map)
#summary(ans2)

####=========================================####
#### compare to GWAS including dominance and epistasis
####=========================================####
#D <- D.mat(marks)
#E <- E.mat(marks)
#ETA <- list( list(Z=Z1, K=K1),list(Z=Z1, K=D),list(Z=Z1, K=E))
#ans3 <- mmer(y=yy, Z=ETA.ADE, W=marks, method="EMMA", 
#              ploidy=4, models=models[1], map=my.map)
#summary(ans2)
}
\keyword{datasets}
