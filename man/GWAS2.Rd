\name{GWAS2}
\alias{GWAS2}
\title{
Genome Wide Association Analysis
}
\description{

Genome Wide Association Analysis using the formula-based multivariate mixed model solver `mmer2`. For details about the generalized least squares (GLS) strategy using the V matrix from the mixed model go to the sommer help page. Please make sure that the marker matrix (W) is sorted correctly and aligned with the phenotype for all individuals. \strong{The function will assume that the observations in the dataset and marker matrix are in the same order}.

}
\usage{

GWAS2(fixed, random, rcov, data, G=NULL, W=NULL, grouping=NULL, method="NR",
      init=NULL,iters=20,tolpar = 1e-06, tolparinv = 1e-06, 
      draw=FALSE, silent=FALSE, constraint=TRUE, EIGEND=FALSE,
      forced=NULL,restrained=NULL,
      
      P3D=TRUE, models="additive", ploidy=2, min.MAF=0.05, 
      gwas.plots=TRUE, map=NULL,manh.col=NULL,
      fdr.level=0.05, n.PC=0)
}
\arguments{
  \item{fixed}{a formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:
  
  \emph{Yield ~ Location} for univariate models
  
  \emph{cbind(Yield,color) ~ Location} for multivariate models.}
  
  \item{random}{a formula specifying the name of the \strong{random effects}, i.e. \emph{random= ~ genotype + year}.
  
  Useful functions can be used to fit heterogeneous variances and others:
  
  \code{\link{at}}(.)   can be used to specify heterogeneous variance, i.e. \emph{random=~at(Location,c("3","4")):ID}
  
  \code{\link{diag}}(.)   can be used to specify a diagonal variance structure, i.e. \emph{random=~diag(Location):ID}
  
  \code{\link{us}}(.)   can be used to specify a unstructured variance, i.e. \emph{random=~us(Location):ID}
  
  \code{\link{and}}(.)   can be used to specify overlay of matrices, i.e. \emph{random=~male + and(female)}
  
  \code{\link{g}}(.)   can be used to specify that a covariance structure for the random effect (not inverse) should be included, i.e. \emph{random=~g(male), G=list(male=Am)}, etc.
  
      \code{\link{grp}}(.)   can be used to specify customized incidence matrices for random effects, i.e. \emph{random=~male + grp(re), grouping=list(re=Zre)}
      
         \strong{For multivariate models} the functions; us(trait) and diag(trait) can be used before the random effects to indicate the covariance structure among traits, i.e.:
  
  \emph{us(trait):Geno} \strong{or} \emph{diag(trait):Geno}
  
  \emph{us(trait):at(Location):Geno} \strong{or} \emph{diag(trait):at(Location):Geno}
  
  \emph{us(trait):us(Location):Geno} \strong{or} \emph{diag(trait):us(Location):Geno}
 
  
  }
  \item{rcov}{a formula specifying the name of the \strong{error term}, i.e. \emph{rcov= ~ units} or \emph{rcov=~at(.):units}.
  
  The functions that can be used to fit heterogeneous residual variances are:
  
  \code{\link{at}}(.)   can be used to specify heterogeneous variance, i.e. \emph{rcov=~at(Location):units}
  
   \strong{For multivariate models} the functions; us(trait) and diag(trait) can be used before the residual term to indicate the covariance structure among traits, i.e.:
  
  \emph{us(trait):units} \strong{or} \emph{diag(trait):units}
  
  \emph{us(trait):at(Location):units} \strong{or} \emph{diag(trait):at(Location):units}
  
  }
  \item{data}{a data frame containing the variables specified in the formulas for response, fixed, and random effects.}
  \item{G}{a list containing the \strong{variance-covariance matrices for the random effects (not inverses)}, i.e. if the user specifies: 
  
  \emph{random=~ g(genotype) + g(year)}
  
  then the G argument should be used as:  
  
  G=list(genotype=M1, year=M2) 
  
  where M1 and M2 are the variance-covariance matrices for random effects 'genotype' and 'year' respectively. Opposite to asreml you provide the original var-covar matrix, not its inverse.}
  
  \item{W}{an incidence matrix for extra fixed effects and only to be used if GWAS is desired and markers will be treated as fixed effects according to Yu et al. (2006) for diploids, and Rosyara et al (2016) for polyploids. Theoretically X and W are both fixed effects, but they are separated to perform GWAS in a model y = Xb + Zu + Wg, allowing the program to recognize the markers from other fixed factors such as environmental factors. This has to be provided as a matrix same than X.
  
  Performs genome-wide association analysis based on the mixed model (Yu et al. 2006):

\deqn{y = X \beta + Z g + W \tau + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The terminology "P3D" (population parameters previously determined) was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the p.adjust function included in the stats package.}

  \item{grouping}{argument to provide customized random effects. This has to be provided as a list where each element of the list is an incidence matrix for a random effect specified in the random formula using the `grp()` function. For example a customized random effect 're' would specified as:
  
  \emph{random=~ grp(re), grouping=list(re=Zre)}
  
  where Zre is the customized incidence matrix. Additional functions such as g() can be used and the variance covariance matrix is provided in the G argument as any other random effect.
  
  }
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. Sommer now focus only in 2 algorithms; Direct-inversion Newton-Raphson \strong{NR} (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2015), and \strong{EMMA} efficient mixed model association (Kang et al. 2008).}
  \item{init}{initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional. It has to be provided as a list, where each list element is one variance component and if multitrait model is pursued each element of the list is a matrix of variance covariance components among traits.}
  \item{iters}{Maximum number of iterations allowed.  Default value is
    15.}
  \item{tolpar}{Convergence criteria.}
  \item{tolparinv}{tolerance parameter for matrix inverse}
    \item{draw}{a TRUE/FALSE value indicating if a plot showing the progress of the variance components should be drawn or not. Default is FALSE}
    \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar as iterations progress. The default is FALSE, which means is not silent and will display the progress bar.}
    \item{constraint}{a TRUE/FALSE value indicating if the function should apply the boundary constraint indicating that variance components that are zero should be restrained.}
    \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this method are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
    \item{forced}{a list of values for variance-covariance components to be used if the user wants to force those values. The format is the same than the init argument.}
   \item{restrained}{a numeric argument specifying which variance-covariance parameters should be restrained. Please see \code{\link{vctable}} and \code{\link{vctable.help}} functions if you want to use it although you shouldn't since is made for the \code{\link{GWAS2}} function.}
   \item{P3D}{when the user performs GWAS, P3D=TRUE means that the variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}
  \item{models}{The model to be used in GWAS. The default is the additive model which applies for diploids and polyploids but the model can be a vector with all possible models, i.e. "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref" models are supported for polyploids based on Rosyara (2016).}
  \item{ploidy}{A numeric value indicating the ploidy level of the organism. The default is 2 which means diploid but higher ploidy levels are supported. This should only be modified if you are performing GWAS in polyploids.}
  \item{min.MAF}{when the user performs GWAS min.MAF is a scalar value between 0-1 indicating what is theminor allele frequency to be allowed for a marker during a GWAS analysis when providing the matrix of markers W. In general is known that results for markers with alleles with MAF < 0.05 are not reliable unless sample size is big enough.}
  \item{gwas.plots}{a TRUE/FALSE statement indicating if the GWAS and qq plot should be drawn or not. The default is TRUE but you may want to avoid it.}
  \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' not neccesarily with same dimensions that markers. The program will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
  \item{manh.col}{a vector with colors desired for the manhattan plot. Default are cadetblue and red alternated.}
  \item{fdr.level}{a level of FDR to calculate and plot the line in the GWAS plot. Default is fdr.level=0.05}
  \item{n.PC}{when the user performs GWAS this refers to the number of principal components to include as fixed effects for Q + K model. Default is 0 (equals K model).}
    
}
\details{

For details about the models enabled please check the help page ot the package \code{\link{sommer}}.

}
\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{var.comp}{a list with the values of the variance-covariance components with one list element for each random effect.}
\item{V.inv}{the inverse of the phenotypic variance matrix V^- = (ZGZ+R)^-1}
\item{u.hat}{a list (one element for each random effect) with a data frame for trait BLUPs.}
\item{Var.u.hat}{a list (one element for each random effect) with the variance-covariance matrix for trait BLUPs.}
\item{PEV.u.hat}{a list (one element for each random effect) with the predicted error variance matrix for trait BLUPs.}
\item{beta.hat}{a data frame for trait BLUEs (fixed effects).}
\item{Var.beta.hat}{a variance-covariance matrix for trait BLUEs}
\item{fish.inv}{inverse of the Fisher's information or average information matrices to obtain variance-covariance of the variance components.}
\item{residuals}{Residual values e = Y - XB}
\item{cond.residuals}{Conditional residual e.cond = Y - (XB + ZU)}
\item{LL}{LogLikelihood}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{X}{ incidence matrix for fixed effects}
\item{dimos}{dimnensions for incidence matrix for random effects}
\item{sigma.scaled}{scaled variance covariance components}
\item{fitted.y}{ Fitted values y.hat=XB+Zu}
\item{fitted.u}{ Fitted values only across random effects u.hat=Zu.1+....+Zu.i}
\item{ZETA}{Original incidence and variance covariance matrices used in the model fit.}
\item{K}{ variance-covariance matrix for random effects. If more than one random effect this will be the diagonal binding of individual K matrices.}
\item{fish.inv}{ If was set to TRUE the Fishers information matrix will be in this slot.}
\item{method}{The method for extimation of variance components specified by the user.}
\item{forced}{a TRUE/FALSE statement indicating if user used the forced argument.}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{W}{The markers matrix used in the GWAS analysis (bad markers discarded).}
\item{W.scores}{A list with two elements; 

\emph{scores}: which are the -log10 p values for each marker

\emph{beta}: which are the GLS estimates; beta= (X'V-X)-X'V-Y. Where beta contains the effect for the intercept (first row) and the effect for the marker (second row). For multiple traits they will be rbinded in sets of two.

}
\item{map}{the map with scores if a map was provided in the GWAS function.}
\item{fdr}{p-value at the false discovery rate (FDR).}

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42:355-360.

}
\author{
Giovanny Covarrubias-Pazaran
}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####

####=========================================####
####=========================================####
#### EXAMPLE 1
#### GWAS in diploids
####=========================================####
####=========================================####

data(CPdata)
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix 
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)
# mix1 <- GWAS2(color~1,
#               random=~g(id), 
#               rcov=~units,
#               G=list(id=A), 
#               W=CPgeno,
#               data=CPpheno)
# summary(mix1)
# 
# ####=========================================####
# ####=========================================####
# #### EXAMPLE 2
# #### GWAS in tetraploids
# ####=========================================####
# ####=========================================####
# 
# data(PolyData)
# genotypes <- PolyData$PGeno
# phenotypes <- PolyData$PPheno
# 
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
# ### run the GWAS model
# ###=========================================####
# my.map <- PolyData$map
# models <- c("additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
# ans2 <- GWAS2(tuber_shape~1, random=~g(Name), models = "additive",
#               G=list(Name=A), W=marks, data=phenotypes2)
# summary(ans2)


}