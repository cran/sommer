\name{BTdata}
\alias{BTdata}
\docType{data}
\title{
Blue Tit Data for a Quantitative Genetic Experiment
}
\description{

a data frame with 828 rows and 7 columns, with variables tarsus length (tarsus) and colour (back) measured on 828 individuals (animal). The mother of each is also recorded (dam) together with the foster nest (fosternest) in which the chicks were reared. The date on which the first egg in each nest hatched (hatchdate) is recorded together with the sex (sex) of the individuals.

}
\usage{data("BTdata")}
\format{
  The format is:
 chr "BTdata"
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples

####=========================================####
####=========================================####
#### EXAMPLE 1
#### simple example
####=========================================####
####=========================================####
data(BTdata)
head(BTdata)
# mix4 <- mmer2(tarsus ~ sex, random = ~ dam + fosternest,
#               data = BTdata) 
# summary(mix4)
#### calculate the ratio and its SE
# pin(mix4, dam.prop ~ V1 / ( V1 + V2 + V3 ) )

#   ####=========================================####
#   ####=========================================####
#   ####=========================================####
#   #### EXAMPLE 2
#   #### more complex multivariate model
#   ####=========================================####
#   ####=========================================####
#   data(BTdata)
#   mix3 <- mmer2(cbind(tarsus, back) ~ sex,
#                 random = ~ dam + fosternest,
#                 data = BTdata, MVM=TRUE)
#   summary(mix3)
#   #### calculate the genetic correlation
#   pin(mix3, gen.cor ~ V2 / sqrt(V1*V3))

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and and \code{\link{mmer2}}}
