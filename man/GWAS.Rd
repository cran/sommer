\name{GWAS}
\alias{GWAS}
\title{
Genome Wide Association Analysis
}
\description{
Genome Wide Association Analysis using the multivariate mixed model solver `mmer`. For details about the generalized least squares (GLS) strategy using the V matrix from the mixed model go to the sommer help page. Please make sure that the marker matrix (W) is sorted correctly and aligned with the phenotypes (Y) for all individuals. \strong{The function will assume that the genotypes in the dataset and marker matrix are in the same order}.
}
\usage{
GWAS(Y, X=NULL, Z=NULL, R=NULL, W=NULL, 
     method="NR", tolpar = 1e-06, tolparinv = 1e-06, 
     draw=FALSE, silent=FALSE, iters=15, 
     init=NULL, che=TRUE, EIGEND=FALSE, forced=NULL, 
     
     P3D=TRUE, models="additive", ploidy=2, min.MAF=0.05, 
     gwas.plots=TRUE, map=NULL,manh.col=NULL,
     fdr.level=0.05, constraint=TRUE, n.PC=0, 
     restrained=NULL)
}
\arguments{
  \item{Y}{a matrix or data frame of response variables}
  \item{X}{an incidence matrix for fixed effects.}
  \item{Z}{incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. \strong{This needs to be provided as a 2-level list structure}. For example:
  
  '
  
  ETA <- list(
  
  A=list(Z=Z1, K=K1),
  
  B=list(Z=Z2, K=K2),
  
  C=list(Z=Z3, K=K3)
  
  ) 
  
  
  '
  
  makes a 2 level list for 3 the random effects A, B and C, stored in a variable we call ETA. The general idea is that each random effect is a list, i.e. A=list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix for the random effect, \strong{if K is not provided is assumed an identity matrix} conferring independence. 
  
  '
  
  PLEASE remember to \strong{use the names Z and K for all random effects} when you provide your matrices, \strong{that's the only way the program distinguishes between a Z or a K matrix}. 
  
  '
  
  To provide extra detail, I'll rephrase it; when moving to situations of more than one random effect, you need to build a list for each random effect, and at the end everything gets joined in a list as well (BGLR type of format). Is called a 2-level list, i.e. A=list(Z=Z1, K=K1) and B=list(Z=Z2, K=K2) refers to 2 random effects and they should be put together in a list:
  
  '
  
  ETA <- list( A=list(Z=Z1, K=K1), B=list(Z=Z1, K=K1) )
  
  '
  
  Now you can fit your model as:
  
  '
  
  mod1 <- mmer(Y=y, Z=ETA) 
  
  '
  
  You can see the examples at the bottom to have a clearer idea how to fit your models. }
 \item{R}{list of R matrices for residual effects.}
 
 \item{W}{an incidence matrix for extra fixed effects and only to be used if GWAS is desired and markers will be treated as fixed effects according to Yu et al. (2006) for diploids, and Rosyara et al (2016) for polyploids. Theoretically X and W are both fixed effects, but they are separated to perform GWAS in a model y = Xb + Zu + Wg, allowing the program to recognize the markers from other fixed factors such as environmental factors. This has to be provided as a matrix same than X.
  
  Performs genome-wide association analysis based on the mixed model (Yu et al. 2006):

\deqn{y = X \beta + Z g + W \tau + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The terminology "P3D" (population parameters previously determined) was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the p.adjust function included in the stats package.}

  \item{method}{this refers to the method or algorithm to be used for estimating variance components. Sommer now focus only in 2 algorithms; Direct-inversion Newton-Raphson \strong{NR} (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2015), and \strong{EMMA} efficient mixed model association (Kang et al. 2008). Average Information has been discontinued because of the instability for very complex models requiring EM-type of algorithms to work.}
  \item{tolpar}{Convergence criteria.}
  \item{tolparinv}{tolerance parameter for matrix inverse}
  \item{draw}{a TRUE/FALSE value indicating if a plot showing the progress of the variance components should be drawn or not. Default is FALSE}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar as iterations progress. The default is FALSE, which means is not silent and will display the progress bar.}
  \item{iters}{Maximum number of iterations allowed.  Default value is
    15.}
  \item{init}{initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional. It has to be provided as a list, where each list element is one variance component and if multitrait model is pursued each element of the list is a matrix of variance covariance components among traits.}
  \item{che}{a TRUE/FALSE value indicating if list structure provided by the user is correct to fix it. The default is TRUE but is turned off to FALSE within the mmer function which would imply a double check.}
 \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this methos are:
  1) can only be applied to one relationship matrix
  2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
  The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
  \item{forced}{a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
  \item{P3D}{when the user performs GWAS, P3D=TRUE means that the variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}
  \item{models}{The model to be used in GWAS. The default is the additive model which applies for diploids and polyploids but the model can be a vector with all possible models, i.e. "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref" models are supported for polyploids based on Rosyara (2016).}
  \item{ploidy}{A numeric value indicating the ploidy level of the organism. The default is 2 which means diploid but higher ploidy levels are supported. This should only be modified if you are performing GWAS in polyploids.}
  \item{min.MAF}{when the user performs GWAS min.MAF is a scalar value between 0-1 indicating what is theminor allele frequency to be allowed for a marker during a GWAS analysis when providing the matrix of markers W. In general is known that results for markers with alleles with MAF < 0.05 are not reliable unless sample size is big enough.}
  \item{gwas.plots}{a TRUE/FALSE statement indicating if the GWAS and qq plot should be drawn or not. The default is TRUE but you may want to avoid it.}
  \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' not neccesarily with same dimensions that markers. The program will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
  \item{manh.col}{a vector with colors desired for the manhattan plot. Default are cadetblue and red alternated.}
  \item{fdr.level}{a level of FDR to calculate and plot the line in the GWAS plot. Default is fdr.level=0.05}
  \item{constraint}{a TRUE/FALSE value indicating if the function should apply the boundary constraint indicating that variance components that are zero should be removed from the analysis and variance components recalculated.}
  \item{n.PC}{when the user performs GWAS this refers to the number of principal components to include as fixed effects for Q + K model. Default is 0 (equals K model).}
  \item{restrained}{a numeric argument specifying which variance-covariance parameters should be restrained. Please see \code{\link{vctable}} and \code{\link{vctable.help}} functions if you want to use it although you shouldn't since is made for the \code{\link{mmer2}} function.}

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

\emph{beta}: which are the GLS estimates of the marker effects, beta= (X'V-X)-X'V-Y.

}
\item{map}{the map with scores if a map was provided in the GWAS function.}
\item{fdr}{p-value at the false discovery rate (FDR).}

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

####=========================================####
#### EXAMPLE 1
#### GWAS in diploids
####=========================================####

data(CPdata)
### look at the data
head(CPpheno)
CPgeno[1:5,1:5]
## fit a model including additive and dominance effects
Y <- CPpheno
Za <- diag(dim(Y)[1])
A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
ETA.A <- list(add=list(Z=Za,K=A))
#ans<- GWAS(Y=Y[,c("color")], Z=ETA.A, W=CPgeno)
#summary(ans)

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

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
# Za <- model.matrix(~Name-1, phenotypes)
# ###=========================================####
# ### run the GWAS model
# ###=========================================####
# head(phenotypes)
# my.map <- PolyData$map
# models <- c("additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref")
# ETA <- list(add=list(Z=Za,K=A))
# ans2 <- GWAS(Y=phenotypes2$tuber_shape,  Z=ETA, 
#              models = "additive", W=marks)
# summary(ans2)


}