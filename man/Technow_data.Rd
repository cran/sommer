\name{Technow_data}
\alias{Technow_data}
\docType{data}
\title{
Genotypic and Phenotypic data from single cross hybrids (Technow et al. (2014))
}
\description{
This dataset contains phenotpic data for 2 traits measured in 1254 single cross hybrids coming from the cross of Flint x Dent heterotic groups. In addition contains the genotipic data (35,478 markers) for each of the 123 Dent lines and 86 Flint lines. The purpose of this data is to demosntrate the prediction of unrealized crosses (9324 unrealized crosses, 1254 evaluated, total 10578 single crosses). We have added the additive relationship matrix (A) but can be easily obtained using the A.mat function on the marker data. Please if using this data for your own research cite Technow et al. (2014) publication (see References).

}
\usage{data("Technow_data")}
\format{
  The format is:
 chr "Technow_data"
}
\source{
This data was extracted from Technow et al. (2014).
}
\references{
If using this data for your own research please cite:

Technow et al. 2014. Genome properties and prospects of genomic predictions of hybrid performance in a Breeding program of maize. Genetics 197:1343-1355.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(Technow_data)

A.flint <- Technow_data$AF # Additive relationship matrix Flint
A.dent <- Technow_data$AD # Additive relationship matrix Dent
M.flint <- Technow_data$MF # Marker matrix Flint
M.dent <- Technow_data$MD # Marker matrix Dent

pheno <- Technow_data$pheno # phenotypes for 1254 single cross hybrids
pheno$hy <- paste(pheno$dent, pheno$flint, sep=":");head(pheno);dim(pheno) 

####=========================================####
#### CREATE A DATA FRAME WITH ALL POSSIBLE HYBRIDS
####=========================================####
# DD <- kronecker(A.dent,A.flint,make.dimnames=TRUE)
# hybs <- data.frame(sca=rownames(DD),yield=NA,matter=NA,gcad=NA, gcaf=NA)
# hybs$yield[match(pheno$hy, hybs$sca)] <- pheno$GY
# hybs$matter[match(pheno$hy, hybs$sca)] <- pheno$GM
# hybs$gcad <- as.factor(gsub(":.*","",hybs$sca))
# hybs$gcaf <- as.factor(gsub(".*:","",hybs$sca))
# head(hybs)
####=========================================####
## RUN THE PREDICTION MODEL
####=========================================####
# y.trn <- hybs
# vv1 <- which(!is.na(hybs$yield))
# vv2 <- sample(vv1, 100)
# y.trn[vv2,"yield"] <- NA
# anss2 <- mmer2(yield~1, random=~g(gcad) + g(gcaf), G=list(gcad=A.dent, gcaf=A.flint), 
#                method="EM", data=y.trn) 
# summary(anss2)
# cor(anss2$fitted.y[vv2], hybs$yield[vv2])

#### try EM algorithm as well if only 2 variance components
#### you can try adding SCA effects by adding SCA=list(Z=Z3, K=DD)
#### to the random effects but will slow down the model a lot


############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 2
#### GWAS in single cross hybrids
####=========================================####
####=========================================####
# data(Technow_data)
# pheno <- Technow_data$pheno # phenotypes for 1254 single cross hybrids
# pheno$hy <- paste(pheno$dent, pheno$flint, sep=":");head(pheno);dim(pheno) 
# pheno$dent <- as.factor(pheno$dent)
# pheno$flint <- as.factor(pheno$flint)
# ###=========================================####
# ### RUN GWAS combining markers
# ###=========================================####
# M.dent <- Technow_data$MD # Marker matrix Dent
# M.flint <- Technow_data$MF # Marker matrix Flint
# Mall <- rbind(M.dent,M.flint); dim(Mall)
# Adf <- A.mat(Mall) # Additive relationship matrix
# # marker matrix for GWAS is the cbind of dent and flint, 
# # parental genomes analyzed by separate but in same framework
# MGF <- cbind(M.dent[pheno$dent,],M.flint[pheno$flint,]); dim(MGF) 
# 
# mox <- mmer2(GY~1, random=~g(dent) + and(g(flint)), 
#              G=list(dent=Adf, flint=Adf), method="EM",
#              W=MGF, data=pheno)
# summary(mox)

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
