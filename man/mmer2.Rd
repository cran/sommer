\name{mmer2}
\alias{mmer2}
\title{Mixed Model Equations in R v2}
\description{
This function solves univariate and multivariate linear mixed models by different methods. It has been implemented to be more user-friendly than \code{\link{mmer}}(which works directly with incidence matrices and variance covariance matrices for each random effect). This version creates the incidence matrices for the user and is more similar to the 'lm' type of functions, with asreml type of specification. Please refer to \code{\link{mmer}} for further details about optimization methods. In general, the methods implemented by this function currently are; "EMMA" efficient mixed model association (Kang et al. 2008), "AI" average information (Gilmour et al. 1995; Lee et al. 2015), "EM" expectation maximization (Searle 1993; Bernardo 2010), and the default "NR" Newton-Raphson (Tunnicliffe 1989). All algorithms are able to deal with multiple responses and multiple variance components; "AI", "EM", "NR", "EMMA"" but "AI" and "NR" are the most efficient for dense structures. Please keep in mind that THIS VERSION NAMED 'mmer2' IS LIMITED compared with \code{\link{mmer}} since the matrices for random effects are built based on a dataframe provided in the 'data' argument, and genetic models sometimes require an incidence matrix with markers (several zero-one columns) or a half diallel matrix,  to count as single random effect. For fitting more general models please use the \code{\link{mmer}} function which works directly with the incidence matrices and variance covariance matrices for each random effect. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices that have been shown to increase prediction accuracy. rrBLUP results can be recreated using the \strong{EMMA} method.

The package has been equiped with several datasets to learn how to use the sommer package; the \code{\link{HDdata}} and \code{\link{FDdata}} datasets will introduce users to fit half and full diallel designs respectively. The \code{\link{h2}} dataset shows how to calculate heritability. The \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets contain data to teach users how to perform genomic prediction in hybrid single crosses which display GCA and SCA effects. The \code{\link{wheatLines}} dataset teaches how to do genomic prediction in single crosses in species displaying only additive effects. The \code{\link{CPdata}} dataset contains data to teach users how to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects.  The same data example is used to show how to include the top GWAS hits as fixed effects in the GBLUP model to increase prediction accuracy, the examples can be found in the \code{\link{hits}} documentation.  The \code{\link{PolyData}} dataset shows how to fit genomic prediction and GWAS analysis in polyploids. The second example in \code{\link{Technow_data}} data shows how to perform GWAS in single cross hybrids. A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. The \code{\link{RICE}} dataset can teach users how to select the best training population using the \code{\link{TP.prep}} function in an applied breeding program, and we show comparison some comparison with other methods. Traces of selection can be obtained using markers with the \code{\link{eigenGWAS}} function. Additional examples for fitting mixed models, such as GWAS and others, can be found in the example section of the \code{\link{mmer}} and \code{\link{mmer2}} functions. Examples of multivariate models can be found in the example #3 of the \code{\link{CPdata}}. An example of how to model autocorrelation structures in the residual components can be found in the \code{\link{yates.oats}} documentation.

PLEASE UPDATE 'sommer' EVERY MONTH USING:

'

install.packages("sommer")

}
\usage{
mmer2(fixed, random, rcov, G=NULL, W=NULL, method="NR", REML=TRUE, 
      MVM=FALSE, iters=20, draw=FALSE, init=NULL, data, family=gaussian, 
      silent=FALSE, constraint=TRUE, sherman=FALSE, EIGEND=FALSE, gss=TRUE,
      forced=NULL, map=NULL, fdr.level=0.05, manh.col=NULL, min.n=FALSE,
      gwas.plots=TRUE,n.cores=1,tolpar = 1e-06,tolparinv = 1e-06)
}
\arguments{
  \item{fixed}{a formula specifying the response variable(s) and the fixed effects, i.e. Yield ~ Location for univariate models, or cbind(color,Yield) ~ Location for multivariate models or univariate in parallel. For running multivariate the 'MVM' argument needs to be set to TRUE.}
  \item{random}{a formula specifying the name of the random effects related to environmental effects, i.e. random= ~ genotype + year. }
    \item{rcov}{a formula specifying the residual structures. The formula has to be specified using the R matrices the user wants to do kronecker. For example, in the spatial designs we want to model the yield of varieties accounting for row and column variation, that could be specified as:
    
    mmer2(yield~1,random=~variety, rcov=~ar1(row):ar1(col), data=your.data)
    
This will create two autocorrelation structures, one for rows and other for columns and do the kroneker to estimate by REML the value of the parameters. The correlation structures available currently are:

"ar1" autocorrelation matrix of order 1

"cs" compound symmetry matrix

"arma" autocorrelation moving average

"id" diagonal matrix

    }
  \item{G}{a list containing the variance-covariance matrices for the random effects, i.e. G=list(genotype=M1, year=M2), where M1 and M2 are the variance-covariance matrices. if not passed is assumed an identity matrix.}
  \item{W}{an incidence matrix for extra fixed effects and only to be used if GWAS is desired and markers will be treated as fixed effects according to Yu et al. (2006) for diploids, and Rosyara et al (2016) for polyploids. Theoretically X and W are both fixed effects, but they are separated to perform GWAS in a model y = Xb + Zu + Wg, allowing the program to recognize the markers from other fixed factors such as environmental factors. This has to be provided as a matrix same than X.
  
  Performs genome-wide association analysis based on the mixed model (Yu et al. 2006):

\deqn{y = X \beta + Z g + W \tau + \varepsilon}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  For unbalanced designs where phenotypes come from different environments, the environment mean can be modeled using the fixed option (e.g., fixed="env" if the column in the pheno data.frame is called "env").  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix.

The terminology "P3D" (population parameters previously determined) was introduced by Zhang et al. (2010).  When P3D=FALSE, this function is equivalent to EMMA with REML (Kang et al. 2008).  When P3D=TRUE, it is equivalent to EMMAX (Kang et al. 2010).  The P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

The dashed line in the Manhattan plots corresponds to an FDR rate of 0.05 and is calculated using the p.adjust function included in the stats package.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. The package currently is supported by 4 algorithms; "EMMA" efficient mixed model association (Kang et al. 2008), "AI" average information (Gilmour et al. 1995; Lee et al. 2015), "EM" expectation maximization (Searle 1993; Bernardo 2010), and "NR" Newton-Raphson (Tunnicliffe 1989). The default method is average information "AI" because of its ability to handle multiple random effects and its greater speed compared to "EM" which can handle multiple random effects but it is much slower.}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{MVM}{a TRUE/FALSE value indicating if the model should be run as multivariate or as parallel univariate models. the default is FALSE which will indicate to run parallel univariate models. The 'n.cores' argument decides how many cores use in such parallelization. MVM=TRUE will run a multivariate model but only for a single random effect.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM or AI algorithms are selected. There is no rule of tumb for the number of iterations. The default value is 50 iterations or EM steps, but could take less or much longer than that. For the AI algorithm usually takes just a few iterations.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is FALSE. COMPUTATION TIME IS SHORTER IF YOU DON'T PLOT SETTING draw=FALSE}
  \item{init}{an vector of initial values for the EM algorithm if this is the method selected. If not provided the program uses a starting values the variance(y)/#random.eff which is usually a good starting value.}
  \item{data}{a data frame containing the variables specified in the formulas for response, fixed, and random effects.}
    \item{family}{a family object to specify the distribution of the response variable. The program will only use the link function to transform the response. For details see \code{\link{family}} help page. The argument would look something like this; family=poisson(), or family=Gamma(), etc. For more sophisticated models please look at lme4 package from Douglas Bates. NOT IMPLEMENTED YET.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar and poems (see \code{\link{poe}} function) while working or should not be displayed. The default is FALSE, which means is not silent and will display the progress bar and a short poem to help the scientist (and me haha) remember that life is more than analyzing data.}
      \item{constraint}{a TRUE/FALSE value indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
      \item{sherman}{a TRUE/FALSE value indicating if Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467) should be used when estimating variance components in order to perform faster when a mixed model with no covariance structure using the average information algorithm is fitted. The default is FALSE since this software was designed for unreplicated data (altough can fit models with replicated data but slower than lme4).}
      \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this methos are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
     \item{gss}{a TRUE/FALSE value indicating if a genomic selection is being fitted just for using certain constraints. When is FALSE (default) the program can make some EM steps to find initial values for variance components when the starting values are to far from the real values causing the likelihood to have a strange behavior and dropping dramatically When TRUE the program does not try EM steps even when far away from the likelihood because in big marker-based models can make the process quite slow.}
    \item{forced}{a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
    \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' not neccesarily with same dimensions that markers. The program will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
    \item{fdr.level}{a level of FDR to calculate and plot the line in the GWAS plot. Default is fdr.level=0.05}
    \item{manh.col}{a vector with colors desired for the manhattan plot. Default are cadetblue and red alternated.}
    \item{min.n}{a TRUE/FALSE statement indicating if the constraint of number of levels of each grouping factor must be < number of observations. This is a constrained usually find in lme4 and has been added and set to TRUE but can be set to FALSE when only one measure per plant is available and the user wants to perform GWAS or genomic selection with limited data.}
    \item{gwas.plots}{a TRUE/FALSE statement indicating if the GWAS and qq plot should be drawn or not. The default is TRUE but you may want to avoid it.}
    \item{n.cores}{number of cores to use when the user passes multiple responses in the model for a faster execution. The default is 1. It relies on forking and hence is not available on Windows unless mc.cores = 1.}
    \item{tolpar}{tolerance parameter for convergence in the multivariate models.}
    \item{tolparinv}{tolerance parameter for matrix inverse in the multivariate models.}
    
}
\details{
The package has been developed to provide R users with code to understand how most common algorithms in mixed model analysis work related to genetics field, but also allowing to perform their real analysis. This package allows the user to calculate the variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects. This program focuses in the mixed model of the form:

'

.............................    Y = Xb + Zu + e   ........ with distributions:

'

Y ~ N[Xb, var(Zu+e)]   ......where;

'

b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, G]  ....... where G is equal to:

'

|K1*sigma2(u1)......................0...........................0.........|
    
|.............0.............K2*sigma2(u2).......................0.........| = G

|......................................................................................|

|.............0....................0.........................Ki*sigma2(ui)...|

'

for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and

'

e ~ N[0, R]  .....................where: I*sigma(e) = R

'

also Var(y) = Var(Zu+e) = ZGZ+R = V which is the phenotypic variance

-------------------------------------------------------------------------------

The functions in the sommer packages allow the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with three random effects:

'

fixed = only intercept...................................b ~ N[b.hat, 0]

random = GCA1 + GCA2 + SCA.................u ~ N[0, G]       

'

then G takes the form:

'

|K1*sigma2(gca1).....................0..........................0.........|
    
|.............0.............K2*sigma2(gca2).....................0.........| = G

|.............0....................0......................K3*sigma2(sca)...|

'

This mixed model would be specified in the \code{\link{mmer}} function as:

'

X1 <- matrix(1,length(y),1)      incidence matrix for intercept only

ETA <- list(

gca1=list(Z=Z1, K=K1), 

gca2=list(Z=Z2, K=K2), 

sca=list(Z=Z3, K=K3)

)     # for 3 random effects

'

where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively created using the \code{\link{model.matrix}} function and K1, K2, K3 are their var-cov matrices. Now the fitted model will be typed as:

'

ans <- mmer(Y=Y, X=X1, Z=ETA)

or

ans <- mmer2(Y~1, random= ~ gca1 + gca2 + sca, G=list(gca1=K1, gca2=K2, sca=K3), data=yourdata)

-------------------------------------------------------------------------------------


PLEASE KEEP IN MIND THAT THIS VERSION NAMED 'mmer2' IS LIMITED compared with \code{\link{mmer}} since the matrices for random effects are built based on a dataframe provided, and genetic models sometimes require a matrix with markers (several columns) to count as single random effect. For fitting more general models please use the \code{\link{mmer}} function which works directly with the incidence matrices and variance covariance matrices for each random effect.

FOR DETAILS ON HOW THE "AI", EM", "NR", AND "EMMA" ALGORITHMS WORK PLEASE REFER TO \code{\link{mmer}}, and methods \code{\link{AI}} , \code{\link{EM}}, \code{\link{EMMA}} AND \code{\link{NR}}

In addition, the package contains a very nice function to plot genetic maps with numeric variable or traits next to the LGs, see the \code{\link{map.plot}} function to see how easy can be done. Other functions contained in the package are:

\code{\link{transp}} function transform a vector of colors in transparent colors.

\code{\link{fdr}} calculates the false discovery rate for a vector of p-values.

\code{\link{A.mat}} is a wrapper of the A.mat function from the rrBLUP package.

\code{\link{score.calc}} is a function that can be used to calculate a -log10 p-value for a vector of BLUEs for marker effects.

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice sommer uses randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using this function which is a wrapper of \code{\link{mmer}}, the core of the sommer package.

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
#### remove them and run the examples
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
}