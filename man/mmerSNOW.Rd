\name{mmerSNOW}
\alias{mmerSNOW}
\title{Mixed Model Equations in R univariate}
\description{
This function is used internally for the mmer function when univariate models are run in parallel. Please refer to the mmer function help page for more detail. 

Finally, feel free to get in touch with me if you have any questions or suggestion at: 

covarrubiasp@wisc.edu

I'll be glad to help or answer any question. We have spent valuable time developing this package. Please cite us in your publication. Type 'citation("sommer")' to know how to cite it.

Please report bugs and provide data to recreate the problem. The only way to improve open-source software is with the scientific community help.

}
\usage{
mmerSNOW(y, X=NULL, Z=NULL, W=NULL, R=NULL, method="AI", REML=TRUE, 
     iters=20, draw=FALSE, init=NULL, n.PC=0, P3D=TRUE,
     models="additive", ploidy=2, min.MAF = 0.05, silent=FALSE, 
     family=NULL, constraint=TRUE, sherman=FALSE, EIGEND=FALSE,
     Fishers=FALSE, gss=TRUE, forced=NULL, full.rank=TRUE, 
     map=NULL,fdr.level=0.05, manh.col=NULL, gwas.plots=TRUE)
}
\arguments{
  \item{y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects related to environmental effects or experimental design. This has to be provided as a matrix, NOT in a list structure.}
  \item{Z}{incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. THIS NEEDS TO BE PROVIDED AS A 2-LEVEL LIST STRUCTURE. For example:
  
  .
  
  ETA <- list(
  
  A=list(Z=Z1, K=K1),
  
  B=list(Z=Z2, K=K2),
  
  C=list(Z=Z3, K=K3)
  
  ) 
  
  
  .
  
  makes a 2 level list for 3 the random effects A, B and C, stored in a variable we call ETA. The general idea is that each random effect is a list, i.e. A=list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix for the random effect, if K is not provided is assumed an identity matrix conferring independence. 
  
  .
  
  PLEASE remember to use the names Z and K FOR ALL RANDOM EFFECTS when you provide your matrices, that's the only way the program distinguishes between a Z or a K matrix. 
  
  .
  
  To provide extra detail, I'll rephrase it; when moving to situations of more than one random effect, you need to build a list for each random effect, and at the end everything gets joined in a list as well (BGLR type of format). Is called a 2-level list, i.e. A=list(Z=Z1, K=K1) and B=list(Z=Z2, K=K2) refers to 2 random effects and they should be put together in a list:
  
  .
  
  ETA <- list( A=list(Z=Z1, K=K1), B=list(Z=Z1, K=K1) )
  
  .
  
  Now you can fit your model as:
  
  .
  
  mod1 <- mmerSNOW(y=y, Z=ETA) 
  
  .
  
  You can see the examples at the bottom to have a clearer idea how to fit your models.
  
  }
  \item{W}{an incidence matrix for extra fixed effects and only to be used if GWAS is desired and markers will be treated as fixed effects according to Yu et al. (2006) for diploids, and Rosyara et al (2016) for polyploids. Theoretically X and W are both fixed effects, but they are separated to perform GWAS in a model y = Xb + Zu + Wg, allowing the program to recognize the markers from other fixed factors such as environmental factors. This has to be provided as a matrix same than X.}
  \item{R}{a matrix for variance-covariance structures for the residuals, i.e. for longitudinal data. if not passed is assumed an identity matrix. THIS PART STILLS IN DEVELOPMENT, NOT FUNCTIONAL YET, it is plan to be implemented in version 1.6.}
  \item{method}{this refers to the method or algorithm to be used for estimating variance components. The package currently is supported by 4 algorithms; "EMMA" efficient mixed model association (Kang et al. 2008), "AI" average information (Gilmour et al. 1995; Lee et al. 2015), "EM" expectation maximization (Searle 1993; Bernardo 2010), and Newton-Raphson "NR" (Tunnicliffe 1989). The default method is average information "AI" because of its ability to handle multiple random effects and its greater speed compared to "EM", "NR" and "EMMA" which can handle multiple random effects but are slower in dense models.}
  \item{REML}{a TRUE/FALSE value indicating if restricted maximum likelihood should be used instead of ML. The default is TRUE.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM or AI algorithms are selected. There is no rule of tumb for the number of iterations. The default value is 50 iterations or EM steps, but usually will take less than that stopping before reaching the maximum number of iterations. For the AI algorithm usually takes just a few iterations.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the log-likelihood should be drawn or not during the optimization process. The default is FALSE. It's been set to FALSE because is less the computation time when the computer doesn't have to draw plots.}
  \item{init}{an vector of initial values for the EM, NR or AI algorithms. If not provided the program uses a starting values the variance(y)/#random.eff which are usually good starting values.}
  \item{n.PC}{when the user performs GWAS this refers to the number of principal components to include as fixed effects for Q + K model. Default is 0 (equals K model).}
  \item{P3D}{when the user performs GWAS, P3D=TRUE means that the variance components are estimated by REML only once, without any markers in the model. When P3D=FALSE, variance components are estimated by REML for each marker separately. The default is the first case.}
    \item{models}{The model to be used in GWAS. The default is the additive model which applies for diploids and polyploids but the model can be a vector with all possible models, i.e. "additive","1-dom-alt","1-dom-ref","2-dom-alt","2-dom-ref" models are supported for polyploids based on Rosyara (2016).}
      \item{ploidy}{A numeric value indicating the ploidy level of the organism. The default is 2 which means diploid but higher ploidy levels are supported. This should only be modified if you are performing GWAS in polyploids.}
  \item{min.MAF}{when the user performs GWAS min.MAF is a scalar value between 0-1 indicating what is theminor allele frequency to be allowed for a marker during a GWAS analysis when providing the matrix of markers W. In general is known that results for markers with alleles with MAF < 0.05 are not reliable unless sample size is big enough.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar and poems (see \code{\link{poe}} function) while working or should not be displayed. The default is FALSE, which means is not silent and will display the progress bar and a short poem to help the scientist (and me haha) remember that life is more than analyzing data.}
  \item{family}{a family object to specify the distribution of the response variable. The program will only use the link function to transform the response. For details see \code{\link{family}} help page. The argument would look something like this; family=poisson(), or family=Gamma(), etc. For more sophisticated models please look at lme4 package from Douglas Bates. NOT IMPLEMENTED YET. Planned for ~ v1.8}
      \item{constraint}{a TRUE/FALSE value indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
      \item{sherman}{a TRUE/FALSE value indicating if Sherman-Morrison-Woodbury formula (Seber, 2003, p. 467) should be used when estimating variance components. This will perform faster when a mixed model with no covariance structures is fitted (only AI algorithm). The default is FALSE since this software was designed for unreplicated data (altough can fit models with replicated data but slower than lme4).}
      \item{EIGEND}{a TRUE/FALSE value indicating if an eigen decomposition for the additive relationship matrix should be performed or not. This is based on Lee (2015). The limitations of this method are:
      1) can only be applied to one relationship matrix
      2) The system needs to be squared and no missing data is allowed (then missing data is imputed with the median).
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM.}
      \item{Fishers}{a TRUE/FALSE value indicating if the program should calculate at the final step and return the inverse of the Fishers Information Matrix.}
      \item{gss}{a TRUE/FALSE value indicating if a genomic selection is being fitted just for using certain constraints. When is FALSE the program can make some EM steps to find initial values for variance components when the starting values are to far from the real values causing the likelihood to have a strange behavior and dropping dramatically When TRUE (default) the program does not try EM steps even when far away from the likelihood because in big marker-based models can make the process quite slow.}
      \item{forced}{a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
      \item{full.rank}{a TRUE/FALSE value indicating if the program should investigate X'X to be full rank to avoid problems when solving the linear system. By default this is TRUE which will display a message in the console to le the user know if the X is full rank or not. and will remove extra columns until full rank condition is met. This could not like some users so it can be desactivated but will return an error anyways at some point due to the fact that X'X is not be invertible. This condition is analyzed once missing data in the response variable 'y' has been removed, not in the original X matrix.}
      \item{map}{a data frame with 2 columns, 'Chrom', and 'Locus' not neccesarily with same dimensions that markers. The program will look for markers in common among the W matrix and the map provided. Although, the association tests are performed for all markers, only the markers in common will be plotted.}
      \item{fdr.level}{a level of FDR to calculate and plot the line in the GWAS plot. Default is fdr.level=0.05}
      \item{manh.col}{a vector with colors desired for the manhattan plot. Default are cadetblue and red alternated.}
      \item{gwas.plots}{a TRUE/FALSE statement indicating if the GWAS and qq plot should be drawn or not. The default is TRUE but you may want to avoid it.}

}
\details{
The package has been developed to provide R users with code to understand how most common algorithms in mixed model analysis work related to genetics field, but also allowing to perform their real analysis. This package allows the user to calculate the variance components for a mixed model with the advantage of specifying the variance-covariance structure of the random effects. This program focuses in the mixed model of the form:

.

.............................    y = Xb + Zu + e   ........ with distributions:

.

y ~ N[Xb, var(Zu+e)]   ......where;

.

b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, G]  ....... where G is equal to:

.

|K1*sigma2(u1)......................0...........................0.........|
    
|.............0.............K2*sigma2(u2).......................0.........| = G

|......................................................................................|

|.............0....................0.........................Ki*sigma2(ui)...|

.

for the i.th random effects, allowing the user to specify the variance covariance structures in the K matrices and

.

e ~ N[0, R]  .....................where: I*sigma(e) = R

.

also Var(y) = Var(Zu+e) = ZGZ+R = V which is the phenotypic variance

-------------------------------------------------------------------------------

The functions in the sommer packages allow the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with three random effects:

.

fixed = only intercept...................................b ~ N[b.hat, 0]

random = GCA1 + GCA2 + SCA.................u ~ N[0, G]       

.

then G takes the form:

.

|K1*sigma2(gca1).....................0..........................0.........|
    
|.............0.............K2*sigma2(gca2).....................0.........| = G

|.............0....................0......................K3*sigma2(sca)...|

.

This mixed model would be specified in the \code{\link{mmerSNOW}} function as:

.

X1 <- matrix(1,length(y),1)      incidence matrix for intercept only

ETA <- list(gca1=list(Z=Z1, K=K1), gca2=list(Z=Z2, K=K2), sca=list(Z=Z3, K=K3))      for 3 random effects

.

where Z1, Z2, Z3 are incidence matrices for GCA1, GCA2, SCA respectively created using the \code{\link{model.matrix}} function and K1, K2, K3 are their var-cov matrices. Now the fitted model will be typed as:

.

ans <- mmerSNOW(y=y, X=X1, Z=ETA)

or

ans <- mmerSNOW2(y~1, random= ~ gca1 + gca2 + sca, G=list(gca1=K1, gca2=K2, sca=K3), data=yourdata)

-------------------------------------------------------------------------------------


FOR DETAILS ON HOW THE "AI", EM" AND "EMMA" ALGORITHMS WORK PLEASE REFER TO \code{\link{AI}} , \code{\link{EM}} AND \code{\link{EMMA}}

In addition, the package contains a very nice function to plot genetic maps with numeric variable or traits next to the LGs, see the \code{\link{map.plot}} function to see how easy can be done. The package contains other functions:

\code{\link{transp}} function transform a vector of colors in transparent colors.

\code{\link{fdr}} calculates the false discovery rate for a vector of p-values.

\code{\link{A.mat}} is a wrapper of the A.mat function from the rrBLUP package.

\code{\link{D.mat}} calculates the dominant relationship matrix.

\code{\link{E.mat}} calculates de epistatic relationship matrix.

\code{\link{score.calc}} is a function that can be used to calculate a -log10 p-value for a vector of BLUEs for marker effects.

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice sommer uses randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using this function which is the core of the sommer package.

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

In addition, we have included a couple of random poems from Latin American writers to help the scientist (an me haha) remember from time to time that life is more than analyzing data. You can always silence this feature by setting the argument silent=TRUE, which will avoid the program to display the poems. If you want to contribute with a poem, phrase or short citation for future versions of sommer, feel free to send it to me to: 

covarrubiasp@wisc.edu

Please share your ideas and code, future generations of scientists can be better if we are not greedy sharing our knowledge. Feel free to use my code for your own software! good luck with your analysis.

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
ETA <- list( list(Z=Z1, K=A.mat(M)))
ans <- mmerSNOW(y=y, Z=ETA, method="EMMA")
summary(ans)

####=========================================####
#### run the same but as GWAS 
#### just add the marker matrix in the argument W
#### markers are fixed effects
####=========================================####

#ans <- mmerSNOW(y=y, Z=ETA, W=M, method="EMMA")
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

#ETA <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=S))
#ans <- mmerSNOW(y=y, X=X1, Z=ETA)
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
#fm3 <- mmerSNOW(y=y, Z=ETA) 
####=========================================####
##### or using mmerSNOW2 would look:
####=========================================####
#fm3 <- mmerSNOW2(fixed=milk ~ 1, random = ~ id + herd, 
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
#CPpheno <- CPdata$pheno
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
#ans.A <- mmerSNOW(y=y.trn, Z=ETA.A)
#cor(ans.A$fitted.y[ww], y[ww], use="pairwise.complete.obs")
####=========================####
#### ADDITIVE-DOMINANT MODEL ####
####=========================####
#ETA.AD <- list(list(Z=Za,K=A),list(Z=Zd,K=D))
#ans.AD <- mmerSNOW(y=y.trn, Z=ETA.AD)
#cor(ans.AD$fitted.y[ww], y[ww], use="pairwise.complete.obs")
### greater accuracy !!!! 4 percent increment!!
### we run 100 iterations, 4 percent increment in general
####===================================####
#### ADDITIVE-DOMINANT-EPISTATIC MODEL ####
####===================================####
#ETA.ADE <- list(list(Z=Za,K=A),list(Z=Zd,K=D),list(Z=Ze,K=E))
#ans.ADE <- mmerSNOW(y=y.trn, Z=ETA.ADE)
#cor(ans.ADE$fitted.y[ww], y[ww], use="pairwise.complete.obs")
#### adding more effects doesn't necessarily increase prediction accuracy!

########## NOTE
## nesting in R is indicated as 
## assume blocks nested in locations
## Loc + Block/Loc
## is the same than
## Loc + Block + Loc:Block

}