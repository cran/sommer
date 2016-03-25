\name{CPdata}
\alias{CPdata}
\docType{data}
\title{
Genotypic and Phenotypic data for a CP population (F1; cross between 2 highly heterozygote individuals; i.e. humans, fruit crops, bredding populations in recurrent selection).
}
\description{
This dataset contains phenotpic data for 363 siblings for an F1 cross. These are averages over 2 environments evaluated for 4 traits; color, yield, fruit average weight, and firmness. The columns in the CPgeno file are the markers whereas the rows are the individuals. The CPpheno data frame contains the measurements for the 363 siblings, and as mentioned before are averages over 2 environments.

}
\usage{data("CPdata")}
\format{
  The format is:
 chr "CPdata"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{
Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. https://cran.rstudio.com/web/packages/sommer/ 
}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
data(CPdata)
CPpheno <- CPdata$pheno
CPgeno <- CPdata$geno
### look at the data
head(CPpheno)
CPgeno[1:5,1:5]
## fit a model including additive and dominance effects
y <- CPpheno$color
Za <- diag(length(y))
Zd <- diag(length(y))
Ze <- diag(length(y))
A <- A.mat(CPgeno) # additive relationship matrix
D <- D.mat(CPgeno) # dominant relationship matrix
E <- E.mat(CPgeno) # epistatic relationship matrix

y.trn <- y # for prediction accuracy
ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
y.trn[ww] <- NA
####================####
#### ADDITIVE MODEL ####
####================####
ETA.A <- list(list(Z=Za,K=A))
#ans.A <- mmer(y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================####
#### ADDITIVE-DOMINANCE MODEL ####
####=========================####
ETA.AD <- list(list(Z=Za,K=A),list(Z=Zd,K=D))
#ans.AD <- mmer(y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
### greater accuracy !!!! 4 percent increment!!
### we run 100 iterations, 4 percent increment in general
####===================================####
#### ADDITIVE-DOMINANCE-EPISTASIS MODEL ####
####===================================####
ETA.ADE <- list(list(Z=Za,K=A),list(Z=Zd,K=D),list(Z=Ze,K=E))
#ans.ADE <- mmer(y=y.trn, Z=ETA.ADE)
#cor(ans.ADE$fitted.y[ww], y[ww], use="pairwise.complete.obs")

#summary(ans.A)
#summary(ans.AD)
#summary(ans.ADE)
#### adding more effects doesn't necessarily increase prediction accuracy!
}
\keyword{datasets}
