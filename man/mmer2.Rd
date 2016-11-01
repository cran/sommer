\name{mmer2}
\alias{mmer2}
\title{Mixed Model Equations in R v2}
\description{

This function solves univariate and multivariate linear mixed models by likelihood methods. It has been \strong{implemented to work with a data frame} which will create the incidence and variance covariance matrices for each random effect. Currently 4 optimization \strong{methods} are supported; \strong{efficient mixed model association (EMMA)} (Kang et al. 2008), \strong{average information (AI)} (Gilmour et al. 1995; Lee et al. 2015), \strong{expectation maximization (EM)} (Searle 1993; Bernardo 2010), and the default \strong{Newton-Raphson (NR)} (Tunnicliffe 1989). All algorithms are able to deal with multiple variance components. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices. Covariance structures for the residual component is currently supported only for univariate Newton Raphson. The package provides flexibility to fit other genetic models such as full and half diallel models as well (see the \code{\link{hdm}} function). The core algorithms in our mixed model solver are the Newton-Raphson and Direct Average Information (Lee et al. 2015) which surpasses in performance the MME-based Average Information (Gilmour et al. 1995) when dense covariance structures are present (GBLUP and GWAS case). When these matrices are sparse (\strong{when no covariance structures}), the MME-based is more effective and we recommend the user to \strong{swith to use the lme4} package. \strong{rrBLUP results can be recreated using the EMMA method. Keep in mind that such method is limited to a single random effect}.

Finally, feel free to get in touch with me if you have any \strong{questions, bug reports and suggestions} at: covarrubiasp@wisc.edu

\strong{For a short tutorial} of how to perform different genetic analysis with sommer please type:

vignette("sommer")

}
\usage{
mmer2(fixed, random, rcov, data, G=NULL, W=NULL, method="NR", REML=TRUE, 
      MVM=FALSE, iters=20, draw=FALSE, init=NULL, family=gaussian, 
      silent=FALSE, constraint=TRUE, sherman=FALSE, EIGEND=FALSE,
      forced=NULL, map=NULL, fdr.level=0.05, manh.col=NULL, min.n=FALSE,
      gwas.plots=TRUE,n.cores=1,tolpar = 1e-06,tolparinv = 1e-06)
}
\arguments{
  \item{fixed}{a formula specifying the \strong{response} variable(s) \strong{and fixed effects}, i.e. Yield ~ Location for \strong{univariate} models, or cbind(color,Yield) ~ Location for \strong{multivariate} models or univariate in parallel. For running multivariate the \strong{'MVM' argument} needs to be set to TRUE.}
  \item{random}{a formula specifying the name of the \strong{random effects}, i.e. random= ~ genotype + year. }
    \item{rcov}{a formula specifying the \strong{residual structures}. For example in the spatial designs we want to model the yield of varieties accounting for row and column variation as an autoregressive process, that could be specified as:
    
    mmer2(yield~1,random=~variety, rcov=~ar1(row):ar1(col), data=your.data)
    
This will create two autocorrelation structures, one for rows and other for columns and do the kroneker to estimate by REML the value of the parameters. Still in development so please be careful and report bugs. The correlation structures available currently are:

\strong{ar1} autocorrelation matrix of order 1

\strong{cs} compound symmetry matrix

\strong{arma} autocorrelation moving average

\strong{id} diagonal matrix

    }
  \item{data}{a data frame containing the variables specified in the formulas for response, fixed, and random effects.}
  \item{G}{a list containing the \strong{variance-covariance matrices for the random effects}, i.e. G=list(genotype=M1, year=M2), where M1 and M2 are the variance-covariance matrices. if not passed is assumed an identity matrix.}
    \item{W}{an incidence matrix for extra fixed effects and \strong{only to be used if GWAS is pursued} to fit a model of the form \deqn{y = X \beta + Z g + W \tau + \varepsilon} Markers will be treated as fixed effects according to Yu et al. (2006) \strong{for diploids}, and Rosyara et al (2016) \strong{for polyploids}. For details about the methods please go to the Details section. The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the p.adjust function included in the stats package.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. The package currently is supported by 4 algorithms; \strong{EMMA} efficient mixed model association (Kang et al. 2008), \strong{AI} average information (Gilmour et al. 1995; Lee et al. 2015), \strong{EM} expectation maximization (Searle 1993; Bernardo 2010), and the default Newton-Raphson \strong{NR} (Tunnicliffe 1989).}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{MVM}{a TRUE/FALSE value \strong{indicating if the model should be run as multivariate} or as parallel univariate models. the default is FALSE which will indicate to run parallel univariate models. The 'n.cores' argument decides how many cores use in such parallelization. MVM=TRUE will run a multivariate model but only for a single random effect.}
  \item{iters}{a scalar value indicating \strong{how many iterations} have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than8 is usually enough. The default value is 20 iterations but usually will take less than that stopping before reaching the maximum number of iterations.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of \strong{updated values for the variance components} and the log-likelihood \strong{should be drawn or not} during the optimization process. The default is FALSE. It's been set to FALSE because is less the computation time when the computer doesn't have to draw plots.}
  \item{init}{an vector of initial values for the EM, NR or AI algorithms. If not provided the program uses a starting values the variance(y)/#random.eff which are usually good starting values.}
 \item{family}{a family object to specify the distribution of the response variable. For details see \code{\link{family}} help page. The argument would look something like this; family=poisson(), or family=Gamma(), etc. For more sophisticated models please look at lme4 package from Douglas Bates. \strong{NOT IMPLEMENTED YET}. Planned for ~ v2.5}
\item{silent}{a TRUE/FALSE value indicating if the function should or shouldn't draw the \strong{progress bar}. The default is FALSE, which means is not silent and will display the progress bar.}
      \item{constraint}{a TRUE/FALSE value indicating if the program should use the \strong{boundary constraint} when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
      \item{sherman}{a TRUE/FALSE value indicating if Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467) should be used when estimating variance components.}
      \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this methos are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
    \item{forced}{a vector (univariate models) or list (multivariate models) of numeric values for variance components including error if the user wants to \strong{force the values of the variance components}. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
      \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' to \strong{create a nice manhattan plot when performing GWAS}. The map does not need to have same dimensions that the marker matrix. The function will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
      \item{fdr.level}{a level of \strong{false discovery rate} to calculate and plot the line in the GWAS plot. Default is fdr.level=0.05. If there's not a single significant marker nothing is returned.}
      \item{manh.col}{a vector with \strong{colors desired for the manhattan plot}. Default are cadetblue and red alternated for each chromosome.}
    \item{min.n}{a TRUE/FALSE statement indicating if the constraint of number of levels of each grouping factor must be < number of observations. This is a constrained usually find in lme4 and has been added and set to TRUE but can be set to FALSE when only one measure per plant is available and the user wants to perform GWAS or genomic selection with limited data.}
    \item{gwas.plots}{a TRUE/FALSE statement indicating if the \strong{manhattan and qq plot should be drawn or not}. The default is TRUE but you may want to avoid it.}
      \item{n.cores}{number of cores to use when the user passes multiple responses in the model for a faster execution. The default is 1. It relies on forking and hence is \strong{not available on Windows} unless mc.cores = 1.}
      \item{tolpar}{tolerance parameter for convergence in the multivariate models.}
      \item{tolparinv}{tolerance parameter for matrix inverse in the multivariate models.}
    
}
\details{

The package has been equiped with several datasets to learn how to use the sommer package; the \code{\link{HDdata}} and \code{\link{FDdata}} datasets will introduce users to fit half and full diallel designs respectively. The \code{\link{h2}} dataset shows how to calculate heritability. The \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets contain data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset teaches how to do genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects.  The same data example is used to show how to include the top GWAS hits as fixed effects in the GBLUP model to increase prediction accuracy, the examples can be found in the \code{\link{hits}} documentation.  The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. The second example in \code{\link{Technow_data}} data shows how to perform GWAS in single cross hybrids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. The \code{\link{RICE}} dataset can teach users how to select the best training population using the \code{\link{TP.prep}} function in an applied breeding program, and we show comparison some comparison with other methods. Traces of selection can be obtained using markers with the \code{\link{eigenGWAS}} function. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions. Examples of multivariate models can be found in the example #3 of the \code{\link{CPdata}}. An example of how to model autocorrelation structures in the residual components can be found in the \code{\link{yates.oats}} documentation. In addition, since v2.3 we have added the \code{\link{nna}} function which does nearest neighbor to create a new data frame adjusted for spatial variation. The \code{\link{ExpDesigns}} data contains several datasets to analyze experimental designs relevant to plant breeding and several detailed examples are available.

=======================================================
=======================================================

\strong{MODELS ENABLED}

The mmer function fits linear mixed models by likelihood methods allowing the used of known covariance structures for random effects. If not provided independence is conferred. This program focuses in the mixed model of the form:

'

.............................    Y = Xb + Zu + e   .....................

..................................... or ................................... 

..... Y = x1b1 + ... + xibi + z1u1 + ... + zjuj + e  ..... 


with distributions:

'

Y ~ N[Xb+Zu, var(Zu+e)]   ......where;

'

b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, G]  ....... where G is equal to:

'

|K1*\eqn{\sigma2}(u1)......................0...........................0.........|
    
|.............0.............K2*\eqn{\sigma2}(u2).......................0.........| = G

|...............................................................................|

|.............0....................0.........................Ki*\eqn{\sigma2}(ui)...|

'

for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and

'

e ~ N[0, R]  .....................where: I*\eqn{\sigma2}(e) = R

'

also Var(Y) = Var(Zu+e) = ZGZ+R = V which is the phenotypic variance. Estimation of variance components are optimized with any of the 4 methods available; NR, AI, EM, EMMA (rrBLUP method).

=======================================================
=======================================================

\strong{GWAS MODELS}

The GWAS models in the sommer package are enabled by using the W argument, which is expected to be a numeric marker matrix. Markers are treated as fixed effects according to the model proposed by Yu et al. (2006) for diploids, and Rosyara et al. (2016) (for polyploids). The matrices X and W are both fixed effects, but they are separated by 2 different arguments to distinguish factors such as environmental and design factors for the argument "X" and markers with "W".
  
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

}
}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. EIGEND: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

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
#### breeding values with 3 variance components
#### for hybrid prediction
####=========================================####
####=========================================####
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K
y <- hybrid2$Yield

####=========================================####
#### Realized IBS relationships for set of parents 1
####=========================================####
K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)     
####=========================================####
#### Realized IBS relationships for set of parents 2
####=========================================####
K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)     
####=========================================####
#### SCA relationship matrix (kronecker)
####=========================================####
S <- kronecker(K1, K2) ; dim(S)   
rownames(S) <- colnames(S) <- levels(hybrid2$SCA)

#head(hybrid2)
#ans <- mmer2(Yield ~ Location, random = ~ GCA1 + GCA2 + SCA, 
#       G=list(GCA1=K1, GCA2=K2, SCA=S),data=hybrid2)
#ans$var.comp
#summary(ans)

####=========================================####
#### please remember THIS FUNCTION IS LIMITED since the matrices for random 
#### effects are built based on a dataframe provided, for more general models 
#### including the genomic analysis please use the 'mmer' function which 
#### works directly with matrices and is more flexible
####=========================================####

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

####=========================================####
#### EXAMPLE 2
#### comparison with lmer, install 'lme4' 
#### and run the code below
####=========================================####

#### lmer cannot use var-cov matrices so we will not 
#### use them in this comparison example

#library(lme4)
#library(sommer)
#data(cornHybrid)
#hybrid2 <- cornHybrid$hybrid
#fm1 <- lmer(Yield ~ Location + (1|GCA1) + (1|GCA2) + (1|SCA),
#        data=hybrid2 )
#out <- mmer2(Yield ~ Location, random = ~ GCA1 + GCA2 + SCA,
#        data=hybrid2)
#summary(fm1)
#summary(out)
#### same BLUPs for GCA1, GCA2, SCA than lme4
#plot(out$cond.residuals, residuals(fm1))
#plot(out$u.hat[[1]], ranef(fm1)$GCA1[,1])
#plot(out$u.hat[[2]], ranef(fm1)$GCA2[,1])
#vv=which(abs(out$u.hat[[3]]) > 0)
#plot(out$u.hat[[3]][vv,], ranef(fm1)$SCA[,1])

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

####=========================================####
#### EXAMPLE 3
#### comparison with lmer, install 'lme4' 
#### and run the code below
####=========================================####

#library(lme4)
#data(sleepstudy)
#fm1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
#summary(fm1)

#out <- mmer2(Reaction ~ Days, random = ~ Subject, data=sleepstudy)
#summary(out)

# plot(out$cond.residuals, residuals(fm1))
# plot(out$u.hat[[1]], ranef(fm1)[[1]][,1])

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

####=========================================####
#### EXAMPLE 4
#### Multivariate model example
####=========================================####

# data(FDdata)
# head(FDdata)
# 
# mix <- mmer2(cbind(stems,pods,seeds)~1, 
#              random=~female+male, MVM=TRUE,data=FDdata)
# summary(mix)
# #### genetic variance covariance
# gvc <- mix$var.comp$female
# #### extract variances (diagonals) and get standard deviations 
# sd.gvc <- as.matrix(sqrt(diag(gvc)))
# #### get possible products sd(Vgi) * sd(Vgi')
# prod.sd <- sd.gvc %*% t(sd.gvc)
# #### genetic correlations cov(gi,gi')/[sd(Vgi) * sd(Vgi')] 
# (gen.cor <- gvc/prod.sd)
# #### pods and seeds have a strong negative genetic covariance (-.79)
# #### more pods, less seeds
}