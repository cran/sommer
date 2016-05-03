\name{bag}
\alias{bag}
\title{Creating a fixed effect matrix with significant GWAS markers}
\description{
This function was designed to create a design matrix with the significant markers found by a GWAS analysis in order to use it in the GBLUP or genomic prediction analysis to increase the predction accuracy. The method is based on Abdollahi Arpanahi et al. (2015) paper were using the top 20 GWAS hits markers and using them as fixed effect, increases the prediction accuracy of the genomic prediction model and denominated bag GBLUP. This phenomena has been explained arguing that the mixed model shrinks too much the effect of markers with big effects and therefore using such markers as fixed effects causes a dramatic increase in the prediction accuracy of a model using them. It has to be noted that for traits of low h2 or with a high number of QTL's with small effects the bagging method doesn't help but actually has a reducing effect in prediction accuracy.
}
\usage{
bag(gwasm, nmar=10, threshold=1, pick=FALSE, method="cluster", only.mark=FALSE)
}
\arguments{
  \item{gwasm}{a GWAS model fitted using mmer}
  \item{nmar}{the number of GWAS hits (markers) to be used for designing the incidence matrix. It finds the markers with maximum significance value and uses them to create the design matrix. The default is the top 10 markers.}
  \item{threshold}{a numeric value indicating the minimum significance value to be used for finding the significant markers. the dedault is 1.}
  \item{pick}{a TRUE/FALSE value indicating if the user prefers to pick the peaks by himself. The default is FALSE leaving the peak selection to one of the two methods available. If set to TRUE R will allow the user to pick the peaks by cliking over the peaks and typing 'Esc' when the user is done selecting the peaks.}
  \item{method}{one of the two methods available; "cluster" performs peak selection by making clusters using k-means (random clusters), whereas "maximum" takes the markers with highest log p.values and select those for designing the model matrix.}
  \item{only.mark}{a TRUE/FALSE statement indicating if the only output should be the marker names or the incidence matrix. By default it returns the incidence matrix but if turned to TRUE will return only the marker names. Useful when want to identify significant markers in different populations.}

}

\value{
If all parameters are correctly indicated the program will return:
\describe{

\item{$X2}{a design matrix with individuals in rows and markers in columns to be used for the GBLUP model based on the GWAS model provided.}

}
}
\references{
Abdollahi Arpanahi R, Morota G, Valente BD, Kranis A, Rosa GJM, Gianola D. 2015. Assessment of bagging GBLUP for whole genome prediction of broiler chicken traits. Journal of Animal Breeding and Genetics 132:218-228.
}
\author{
Giovanny Covarrubias-Pazaran
}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with a single '#' mark, 
#### remove them and run the examples
####=========================================####
data(CPdata)
CPpheno <- CPdata$pheno
CPgeno <- CPdata$geno

####=========================================####
#### convert markers to numeric format
####=========================================####
## fit a model including additive and dominance effects
y <- CPpheno$color
Za <- diag(length(y))
A <- A.mat(CPgeno) # additive relationship matrix

####=========================================####
#### identify major genes and create the bagging matrix
####=========================================####

ETA.A <- list(list(Z=Za,K=A))
#ans.GWAS <- mmer(y=y, Z=ETA.A, W=CPgeno)
#summary(ans.GWAS)

####=========================================####
#### run the bag function to design the matrix
#### for top GWAS hits
####=========================================####

#X1 <- bag(ans.GWAS);head(X1); dim(X1)

####=========================================####
#### compare prediction accuracies between
#### GBLUP and bag GBLUP 
####=========================================####
set.seed(1234)
y.trn <- y # for prediction accuracy
ww <- sample(c(1:dim(Za)[1]),72) # delete data for one fifth of the population
y.trn[ww] <- NA

ETA.A <- list(list(Z=Za,K=A))
#ans.A <- mmer(y=y.trn, Z=ETA.A) # GBLUP
#ans.AF <- mmer(y=y.trn, X=X1, Z=ETA.A) # bagging-GBLUP
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs") # GBLUP
#cor(ans.AF$fitted.y[ww], y[ww], use="pairwise.complete.obs") # bagging-GBLUP
#### 11 percent increase in prediction accuracy
}