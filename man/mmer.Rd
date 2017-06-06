\name{mmer}
\alias{mmer}
\title{Mixed Model Equations in R}
\description{
This function solves univariate and multivariate linear mixed models by likelihood methods. It has been \strong{implemented to work directly with incidence and variance covariance matrices} for each random effect. Current optimization \strong{methods} are; \strong{efficient mixed model association (EMMA)} (Kang et al. 2008), \strong{Direct-Inversion Average Information (AI)} (Lee et al. 2015), \strong{MME-based expectation maximization (EM)} (Searle 1993; Bernardo 2010), \strong{MME-based Average Information (AI)}(Gilmour et al. 1995). and the default \strong{Direct-Inversion Newton-Raphson (NR)} (Tunnicliffe 1989). All algorithms are able to deal with multiple variance components. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices. The package provides flexibility to fit other genetic models such as full and half diallel models as well (see how to use the \code{\link{and}} functionality), and heterogeneous variances (\code{\link{at}} function). The core algorithms in sommer based in \strong{Direct-Inversion surpasses in performance the MME-based only when dense covariance structures are present} (GBLUP and GWAS case). \strong{With sparse matrices} (meaning when dense covariance structures don't exist), the MME-based algorithms are faster and we recommend to shift to the use of method="EM" or method="AI" combined with DI="FALSE" to change from direct-inversion to MME-based algorithms. \strong{rrBLUP results can be recreated using the EMMA method. Keep in mind that such method is limited to a single random effect}.

Finally, feel free to get in touch with me if you have any \strong{questions, bug reports and suggestions} at: covarrubiasp@wisc.edu

\strong{For a short tutorial} of how to perform different genetic analysis with sommer please type:

vignette("sommer")

}
\usage{
mmer(Y, X=NULL, Z=NULL, R=NULL, W=NULL, method="NR", REML=TRUE, DI=TRUE,
     MVM=FALSE, iters=20, draw=FALSE, init=NULL, n.PC=0, P3D=TRUE,
     models="additive", ploidy=2, min.MAF = 0.05, silent=FALSE, 
     family=NULL, constraint=TRUE, sherman=FALSE, EIGEND=FALSE,
     forced=NULL, map=NULL,fdr.level=0.05, manh.col=NULL, 
     gwas.plots=TRUE, n.cores=1,tolpar = 1e-06,tolparinv = 1e-06,
     che=TRUE, IMP=TRUE)
}
\arguments{
  \item{Y}{a numeric \strong{vector or matrix/dataframe} for the response variable(s) (it can run multivariate and also runs univariate in parallel by using the 'n.cores' argument). The \strong{multivariate version is only enabled when the argument 'MVM' is set to TRUE} (not the default).}
  \item{X}{an incidence matrix or a list of incidence matrices (in the case of parallel models) for \strong{fixed effects} related to environmental effects or experimental design. These can be constructed using the \code{\link{model.matrix}} function.}
  \item{Z}{incidence matrices and var-cov matrices for \strong{random effects}. This works for ONE OR MORE random effects. \strong{This needs to be provided as a 2-level list structure}. For example:
  
  '
  
  ETA <- list(
  
  A=list(Z=Z1, K=K1),
  
  B=list(Z=Z2, K=K2),
  
  C=list(Z=Z3, K=K3)
  
  ) 
  
  '
  
  mod1 <- mmer(y=y, Z=ETA) 
  
  '
  
  makes a 2 level list for modeling 3 the random effects A, B and C, stored in a variable we named ETA (can have any name). The general idea is that each random effect is a list, i.e. A=list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix for the random effect, \strong{if K is not provided is assumed an identity matrix} confering independence. 
  
  '
  
  PLEASE remember to \strong{use only the names Z and K for all random effects} when you provide your matrices, \strong{that's the only way the program distinguishes between a Z (incidence) or a K (covariance) matrix}. You can see the examples at the bottom to have a clearer idea how to fit your models.
  
  }
  \item{R}{a list with R matrices used when the user wants heterogeneous residual variances, i.e. if user wants different residual variance components for each environment in a multienvironment analysis.}
  \item{W}{an incidence matrix for extra fixed effects and \strong{only to be used if GWAS is pursued} to fit a model of the form \deqn{y = X \beta + Z g + W \tau + \varepsilon}. Markers will be treated as fixed effects according to Yu et al. (2006) \strong{for diploids}, and Rosyara et al (2016) \strong{for polyploids}. For details about the methods please go to the Details section. The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the p.adjust function included in the stats package.}

  \item{method}{this refers to the method or algorithm to be used for estimating variance components. The package currently is supported by 4 algorithms; \strong{EMMA} efficient mixed model association (Kang et al. 2008), \strong{AI} average information (Gilmour et al. 1995; Lee et al. 2015), \strong{EM} expectation maximization (Searle 1993; Bernardo 2010), and the default Newton-Raphson \strong{NR} (Tunnicliffe 1989).}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{DI}{a TRUE/FALSE value indicating if the method (i.e. AI, NR, etc.) to estimate variance components should use direct inversion (DI=TRUE) or MME-based (DI=FALSE) techniques if available for such method.}
  \item{MVM}{a TRUE/FALSE value \strong{indicating if the model should be run as multivariate} or as parallel univariate models. the default is FALSE which will indicate to run parallel univariate models. The 'n.cores' argument decides how many cores use in such parallelization. MVM=TRUE will run a multivariate model but only for a single random effect.}
  \item{iters}{a scalar value indicating \strong{how many iterations} have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than8 is usually enough. The default value is 20 iterations but usually will take less than that stopping before reaching the maximum number of iterations.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of \strong{updated values for the variance components} and the log-likelihood \strong{should be drawn or not} during the optimization process. The default is FALSE. It's been set to FALSE because is less the computation time when the computer doesn't have to draw plots.}
  \item{init}{an vector of initial values for the EM, NR or AI algorithms. If not provided the program uses a starting values the variance(y)/#random.eff which are usually good starting values.}
  \item{n.PC}{when the user performs GWAS this refers to the number of principal components to include as fixed effects for Q + K model. Default is 0 (equals K model). See Details section.}
  \item{P3D}{when the user performs GWAS, P3D=TRUE means that the variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case. See Details section.}
    \item{models}{\strong{The model to be used in GWAS based on ploidy}. The default is the additive model which applies for diploids and polyploids but the model can be a vector with all possible models, i.e. "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref" models are supported for polyploids based on Rosyara (2016). See Details section.}
      \item{ploidy}{A numeric value indicating the ploidy level of the organism. The default is 2 which means diploid but higher ploidy levels are supported. This \strong{is only relevant if you are performing GWAS in polyploids}.}
  \item{min.MAF}{when the user performs GWAS min.MAF is a scalar value between 0-1 indicating what is theminor allele frequency to be allowed for a marker provided in the argument "W". In general is known that results for markers with alleles with MAF < 0.05 are not reliable unless sample size is big enough.}
  \item{silent}{a TRUE/FALSE value indicating if the function should or shouldn't draw the \strong{progress bar}. The default is FALSE, which means is not silent and will display the progress bar.}
  \item{family}{a family object to specify the distribution of the response variable. For details see \code{\link{family}} help page. The argument would look something like this; family=poisson(), or family=Gamma(), etc. For more sophisticated models please look at lme4 package from Douglas Bates. \strong{NOT IMPLEMENTED YET}. Planned for ~ v2.5}
      \item{constraint}{a TRUE/FALSE value indicating if the program should use the \strong{boundary constraint} when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
      \item{sherman}{a TRUE/FALSE value indicating if Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467) should be used when estimating variance components.}
      \item{EIGEND}{a TRUE/FALSE value indicating if an \strong{eigen decomposition for the additive relationship matrix} should be performed or not. This is based on Lee (2015). The limitations of this method are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
      \item{forced}{a vector (univariate models) or list (multivariate models) of numeric values for variance components including error if the user wants to \strong{force the values of the variance components}. On the meantime only works for forcing all of them and not a subset. The default is NULL, meaning that variance components will be estimated by REML/ML.}
      \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' to \strong{create a nice manhattan plot when performing GWAS}. The map does not need to have same dimensions than the marker matrix. The function will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
      \item{fdr.level}{a level of \strong{false discovery rate} to calculate and plot the line \strong{in} the \strong{GWAS} plot. Default is fdr.level=0.05. If there's not a single significant marker nothing is returned.}
      \item{manh.col}{a vector with \strong{colors desired for the manhattan plot}. Default are cadetblue and red alternated for each chromosome.}
      \item{gwas.plots}{a TRUE/FALSE statement indicating if the \strong{manhattan and qq plot should be drawn or not}. The default is TRUE but you may want to avoid it.}
      \item{n.cores}{number of cores to use when the user passes multiple responses in the model for a faster execution of univariate models. The default is 1. It relies on forking and hence is \strong{not available on Windows} unless mc.cores = 1.}
      \item{tolpar}{tolerance parameter for convergence in the models.}
      \item{tolparinv}{tolerance parameter for matrix inversion in the models.}
     \item{che}{a TRUE/FALSE value statement indicating if the program should spend time \strong{checking the matrices} of random effects provided by the user \strong{for possible errors} in dimensions.}
      \item{IMP}{a TRUE/FALSE statement if the function should impute the Y matrix/dataframe with the median value or should get rid of missing values for the estimation of variance components. \strong{Only for multivariate mixed models}.}

}
\details{

The package has been equiped with several datasets to learn how to use the sommer package; the \code{\link{HDdata}} and \code{\link{FDdata}} datasets will introduce users to fit half and full diallel designs respectively. The \code{\link{h2}} dataset shows how to calculate heritability. The \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets contain data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset teaches how to do genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects.  The same data example is used to show how to include the top GWAS hits as fixed effects in the GBLUP model to increase prediction accuracy, the examples can be found in the \code{\link{hits}} documentation.  The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. The second example in \code{\link{Technow_data}} data shows how to perform GWAS in single cross hybrids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. The \code{\link{RICE}} dataset can teach users how to select the best training population using the \code{\link{TP.prep}} function in an applied breeding program, and we show comparison some comparison with other methods. Traces of selection can be obtained using markers with the \code{\link{eigenGWAS}} function. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions. Examples of multivariate models can be found in the example #3 of the \code{\link{CPdata}}. In addition, since v2.3 we have added the \code{\link{nna}} function which does nearest neighbor to create a new data frame adjusted for spatial variation. The \code{\link{ExpDesigns}} data contains several datasets to analyze experimental designs relevant to plant breeding and several detailed examples are available. Useful functions for analyzing such designs are included such as the \code{\link{blocker}} and \code{\link{fill.design}}.

=======================================================
=======================================================

\strong{MODELS ENABLED}

The mmer function fits linear mixed models by likelihood methods allowing the used of known covariance structures for random effects. If not provided independence is conferred. This program focuses in the mixed model of the form:

'

\deqn{Y = X \beta + Z u + \epsilon}

'

with distributions:

'

\deqn{Y ~ MVN ( X\beta+Zu, var(Z u + \epsilon) )}

where;

'

\deqn{\beta  ~ N (\beta, 0)} 

\deqn{u  ~ N (0, G)} 

where G is equal to:

'

\tabular{rrrrr}{
  K1*\eqn{\sigma2}(u1)  \tab \tab 0 \tab \tab 0\cr
  0 \tab \tab K2*\eqn{\sigma2}(u2) \tab \tab 0\cr
  ... \tab \tab ... \tab \tab ...\cr
  0 \tab \tab 0 \tab \tab Ki*\eqn{\sigma2}(ui)\cr
}


'

for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and

'

\deqn{\epsilon  ~ N (0, R)}

where: \eqn{R = I * \sigma2 \epsilon}

'

also Var(Y) = Var(Zu+e) = ZGZ+R = V which is the phenotypic variance. Estimation of variance components are optimized with any of the 4 methods available; NR, AI, EM, EMMA (rrBLUP method).

=======================================================
=======================================================

\strong{GWAS MODELS}

The GWAS models in the sommer package are enabled by using the W argument, which is expected to be a numeric marker matrix. Markers are treated as fixed effects according to the model proposed by Yu et al. (2006) for diploids, and Rosyara et al. (2016) (for polyploids). The matrices X and W are both fixed effects, but they are separated by 2 different arguments to distinguish factors such as environmental and design factors for the argument "X" and markers with "W". Multivariate GWAS models are based on Covarrubias-Pazaran et al. (2017) where covariance among traits are exploited as well.
  
  The genome-wide association analysis is based on the mixed model:

\deqn{y = X \beta + Z g + W \tau + e}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The argument "P3D" was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

For extra details about the methods please read the canonical papers listed in the References section.

=======================================================
=======================================================

Additional functions for genetic analysis have been included such as False Discovery Rate calculation (\code{\link{fdr}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}), phasing F1 or CP populations (\code{\link{phase.F1}}). We have also included the famous \code{\link{transp}} function to get transparent colors. 

}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{

\item{$var.comp}{ a vector with the values of the variance components estimated}
\item{$V.inv}{ a matrix with the inverse of the phenotypic variance V = ZGZ+R, V^-1}
\item{$u.hat}{ a vector with BLUPs for random effects}
\item{$Var.u.hat}{ a vector with variances for BLUPs}
\item{$PEV.u.hat}{ a vector with predicted error variance for BLUPs}
\item{$beta.hat}{ a vector for BLUEs of fixed effects}
\item{$Var.beta.hat}{ a vector with variances for BLUEs}
\item{$LL}{ LogLikelihood}
\item{$AIC}{ Akaike information criterion}
\item{$BIC}{ Bayesian information criterion}
\item{$X}{ incidence matrix for fixed effects}
\item{$fitted.y}{ Fitted values y.hat=XB+Zu}
\item{$fitted.u}{ Fitted values only across random effects u.hat=Zu.1+....+Zu.i}
\item{$residuals}{ Residual values e = y - XB or y - y.hat.fixed}
\item{$cond.residuals}{ Conditional residual values e = y - (XB + Zu) or y - y.hat}
\item{$fitted.y.good}{ Fitted values y.hat=XB+Zu only for data that had no missing data originally. Only used for my checks.}
\item{$Z}{ incidence matrix for random effects. If more than one random effect this will be the column binding of individual Z matrices.}
\item{$K}{ variance-covariance matrix for random effects. If more than one random effect this will be the diagonal binding of individual K matrices.}
\item{$fish.inv}{ If was set to TRUE the Fishers information matrix will be in this slot.}
\item{$method}{ The method for extimation of variance components specified by the user.}
\item{$maxim}{ Maximization used. An argument for the program to know if REML or ML was used. If TRUE means that REML was used instead of ML.}
\item{$score}{ the -log10(p-value) for each marker if a GWAS model is fitted by specifying the W parameter in the model.}
\item{$map}{ if GWAS is performed and a map is provided the program will return a new map of the markers in common among the map and the W matrix and the -log10(p.values) for such marker tests.}

Please share your ideas and code, future generations of scientists can be better trained if we are not greedy sharing our knowledge. Feel free to use my code for your own software! good luck with your analysis.

}
}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744 

Covarrubias-Pazaran G et al. Validating multivariate genomic selection and genome wide association in American cranberry. Submitted to the Plant Genome (2016).

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
#### breeding values with 1 variance component
####=========================================####
####=========================================####

####=========================================####
#### simulate genotypic data
#### random population of 200 lines with 1000 markers
####=========================================####
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
####=========================================####
#### simulate phenotypes
####=========================================####
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
M <- t(M)
####=========================================####
#### fit the model
####=========================================####
Z1 <- diag(length(y))
ETA <- list( add=list(Z=Z1, K=A.mat(M)) )
ans <- mmer(Y=y, Z=ETA)
summary(ans)

####=========================================####
#### run the same but as GWAS 
#### just add the marker matrix in the argument W
#### markers are fixed effects
####=========================================####

#ans <- mmer(Y=y, Z=ETA, W=M)
#summary(ans)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 2
#### breeding values with 3 variance components
#### Hybrid prediction
####=========================================####
####=========================================####
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K
y <- hybrid2$Yield
X1 <- model.matrix(~ Location, data = hybrid2);dim(X1)
Z1 <- model.matrix(~ GCA1 -1, data = hybrid2);dim(Z1)
Z2 <- model.matrix(~ GCA2 -1, data = hybrid2);dim(Z2)
Z3 <- model.matrix(~ SCA -1, data = hybrid2);dim(Z3)

colnames(Z1) <- levels(hybrid2$GCA1)
colnames(Z2) <- levels(hybrid2$GCA2)
colnames(Z3) <- levels(hybrid2$SCA)
####=========================================####
#### Realized IBS relationships for set of parents 1
####=========================================####
#K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1) 
####=========================================####
#### Realized IBS relationships for set of parents 2
####=========================================####
#K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)
####=========================================####
#### Realized IBS relationships for cross 
#### (as the Kronecker product of K1 and K2)
####=========================================####
#S <- kronecker(K1, K2) ; dim(S)   
#rownames(S) <- colnames(S) <- levels(hybrid2$SCA)

#ETA <- list(GCA1=list(Z=Z1, K=K1), 
#            GCA2=list(Z=Z2, K=K2), 
#            SCA=list(Z=Z3, K=S)
#            )
#ans <- mmer(Y=y, X=X1, Z=ETA)
#ans$var.comp
#summary(ans)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 3
#### use data from example 2
#### COMPARE WITH MCMCglmm
####=========================================####
####=========================================####

####=========================================####
#### the same model run in MCMCglmm:
####=========================================####
#library(MCMCglmm)
# pro <- list(GCA1 = as(solve(K1), "sparseMatrix"), GCA2 = as(solve(K2),
#      + "sparseMatrix"), SCA = as(solve(S), "sparseMatrix") )
#system.time(mox <- MCMCglmm(Yield ~ Location, random = ~ GCA1 + GCA2 + SCA,
#      + data = hybrid2, verbose = T, ginverse=pro))
## Takes 7:13 minutes in MCMCglmm, in sommer only takes 7 seconds

####=========================================####
#### it is also possible to do GWAS for hybrids, separatting 
#### and accounting for effects of GCA1, GCA2, SCA
####=========================================####

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 4
#### COMPARE WITH cpgen
####=========================================####
####=========================================####

#Z_list = list(Z1,Z2,Z3)
#G_list = list(solve(K1), solve(K2), solve(S))
#fit <- clmm(y = y, Z = Z_list, ginverse=G_list, niter=15000, burnin=5000)
####=========================================####
#### inspect results and notice that variance 
#### components were NOT estimated correctly!!
#### also takes longer and no user-friendly 
####=========================================####
#str(fit)

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 5
#### COMPARE WITH pedigreemm example
####=========================================####
####=========================================####

#library(pedigreemm)
#A <- as.matrix(getA(pedCowsR))
#y <- milk$milk
#Z1 <- model.matrix(~id-1, data=milk); dim(Z1)
#vv <- match(unique(milk$id), gsub("id","",colnames(Z1)))
#K1<- A[vv,vv]; dim(K1) 
#Z2 <- model.matrix(~as.factor(herd)-1, data=milk); dim(Z2)
#ETA<- list(list(Z=Z1, K=K1),list(Z=Z2))
#fm3 <- mmer(Y=y, Z=ETA) 
####=========================================####
##### or using mmer2 would look:
####=========================================####
#fm3 <- mmer2(fixed=milk ~ 1, random = ~ id + herd, 
#             G=list(id=K1), data=milk)
#summary(fm3)
####=========================================####
#### Try pedigreemm but takes longer, 
#### is an extension of lme4
####=========================================####
#fm2 <- pedigreemm(milk ~ (1 | id) + (1 | herd),data = milk, pedigree = list(id= pedCowsR))
#plot(fm3$u.hat[[1]], ranef(fm2)$id[,1])
#plot(fm3$u.hat[[2]], ranef(fm2)$herd[,1])
####=========================================####
#### a big data frame with 3397 rows and 1359 animals analyzed
#### pedigreemm takes 4 min, sommer takes 1 minute
####=========================================####

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 6
#### PREDICTING SPECIFIC PERFORMANCE 
#### within biparental population    
####=========================================####
####=========================================####

#data(CPdata)
#CPpheno <- CPdata$pheno[,-c(1:4)]
#CPgeno <- CPdata$geno
## look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
####=========================================####
#### fit a model including additive and dominance effects
####=========================================####
#y <- CPpheno$color
#Za <- diag(length(y))
#Zd <- diag(length(y))
#A <- A.mat(CPgeno)
#D <- D.mat(CPgeno)

#y.trn <- y # for prediction accuracy
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#y.trn[ww] <- NA

####================####
#### ADDITIVE MODEL ####
####================####
#ETA.A <- list(list(Z=Za,K=A))
#ans.A <- mmer(Y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================####
#### ADDITIVE-DOMINANT MODEL ####
####=========================####
#ETA.AD <- list(list(Z=Za,K=A),list(Z=Zd,K=D))
#ans.AD <- mmer(Y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
### greater accuracy !!!! 4 percent increment!!
### we run 100 iterations, 4 percent increment in general
####===================================####
#### ADDITIVE-DOMINANT-EPISTATIC MODEL ####
####===================================####
#ETA.ADE <- list(list(Z=Za,K=A),list(Z=Zd,K=D),list(Z=Ze,K=E))
#ans.ADE <- mmer(Y=y.trn, Z=ETA.ADE)
#cor(ans.ADE$fitted.y[ww], y[ww], use="pairwise.complete.obs")
#### adding more effects doesn't necessarily increase prediction accuracy!

########## NOTE
## nesting in R is indicated as 
## assume blocks nested in locations
## Loc + Block/Loc
## is the same than
## Loc + Block + Loc:Block

############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 7
#### MULTIVARIATE MODEL
####=========================================####
####=========================================####
#data(CPdata)
#CPpheno <- CPdata$pheno[,-c(1:4)]
#CPgeno <- CPdata$geno
#### look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
#### fit a model including additive and dominance effects
#Y <- CPpheno
#Za <- diag(dim(Y)[1])
#A <- A.mat(CPgeno) # additive relationship matrix
####================####
#### ADDITIVE MODEL ####
####================####
#ETA.A <- list(add=list(Z=Za,K=A))
#### MAKE SURE YOU SET 'MVM=TRUE'
#ans.A <- mmer(Y=Y, Z=ETA.A, MVM=TRUE)
#summary(ans.A)

#### genetic variance covariance
#gvc <- ans.A$var.comp$Vu
#### extract variances (diagonals) and get standard deviations
#sd.gvc <- as.matrix(sqrt(diag(gvc))) 
#### get possible products sd(Vgi) * sd(Vgi')
#prod.sd <- sd.gvc \%*\% t(sd.gvc) ###remove the \
#### genetic correlations cov(gi,gi')/[sd(Vgi) * sd(Vgi')]
#(gen.cor <- gvc/prod.sd)
#### heritabilities
#(h2 <- diag(gvc) / diag(cov(Y, use = "complete.obs")))


############################################
############################################
############################################
############################################
############################################
############################################

####=========================================####
####=========================================####
#### EXAMPLE 8
#### Genomic prediction
#### Univariate vs Multivariate models
####=========================================####
####=========================================####
data(CPdata)
CPpheno <- CPdata$pheno[,-c(1:4)]
CPgeno <- CPdata$geno
#### look at the data
#head(CPpheno)
#CPgeno[1:5,1:5]
#### fit a model including additive and dominance effects
#Za <- diag(dim(CPpheno)[1])
#A <- A.mat(CPgeno) # additive relationship matrix

#CPpheno2 <- CPpheno
#ww <- sample(c(1:dim(Za)[1]),72) # delete data for 1/5 of the population
#CPpheno2[ww,"Yield"] <- NA
####==================####
#### univariate model ####
####==================####
#ETA.A <- list(add=list(Z=Za,K=A))
#ans.A <- mmer(Y=CPpheno2$Yield, Z=ETA.A, method="NR")
#ans.A$var.comp
#cor(ans.A$fitted.y[ww], CPpheno[ww,"Yield"], use="pairwise.complete.obs")
####====================####
#### multivariate model ####
####====================####
#### estimate var.comp for subset
#ETA.B <- list(add=list(Z=diag(dim(A[-ww,-ww])[1]),K=A[-ww,-ww]))
#ans.B <- mmer(Y=CPpheno2[-ww,c(2,1,3)], Z=ETA.B,MVM=TRUE)
#ans.B$var.comp
#### now force the variance components for getting the plants with no data
#ETA.C <- list(add=list(Z=Za,K=A))
#ans.C <- mmer(Y=CPpheno2[,c(2,1,3)], Z=ETA.C, forced=ans.B$var.comp,MVM=TRUE)
#cor(ans.C$u.hat$add[ww,1], CPpheno[ww,"Yield"], use="pairwise.complete.obs")

}