\name{mmer2}
\alias{mmer2}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations in \strong{R}

}
\description{

Fits a multivariate or univariate linear mixed model by likelihood methods (REML). It has been \strong{implemented to work with a data frame in a formula-based fashion} which will create the incidence and variance covariance matrices for the user and call the \code{\link{mmer}} function in the background. The optimization methods are; Direct-Inversion Newton-Raphson (NR) (Tunnicliffe 1989; Lee et al. 2015; Maier et al. 2015), and Efficient Mixed Model Association (EMMA) (Kang et al. 2008). These algorithms are \strong{intended to be used for problems of the type p > n}. For problems with sparse covariance structures, or problems of the type n > p, the MME-based algorithms are faster and we recommend to shift to the use of MME-based software (i.e. lme4, breedR, or asreml-R).

The package also provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices to model covariances among genotypes for plant and animal breeding problems. 

\strong{For a short tutorial} of how to perform different genetic analysis with sommer please \strong{type}:

\strong{vignette("sommer.start")} or \strong{vignette("sommer")}

}
\usage{

mmer2(fixed, random, rcov, data, G=NULL, grouping=NULL, method="NR",
      init=NULL,iters=20,tolpar = 1e-06, tolparinv = 1e-06, 
      draw=FALSE, silent=FALSE, constraint=TRUE, EIGEND=FALSE,
      forced=NULL, complete=TRUE, restrained=NULL,
      na.method.X="exclude", na.method.Y="exclude", REML=TRUE)
}
\arguments{
  \item{fixed}{a formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:
  
  \emph{Yield ~ Location} for univariate models
  
  \emph{cbind(Yield,color) ~ Location} for multivariate models.}
  
  \item{random}{a formula specifying the name of the \strong{random effects}, i.e. \emph{random= ~ genotype + year}.
  
  Useful functions can be used to fit heterogeneous variances and other special models:
  
  \code{\link{at}}(.)   can be used to specify heterogeneous variance, i.e. \emph{random=~at(Location,c("A","B")):ID}
  
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
   The default is FALSE to avoid the user get into trouble but experimented users can take advantage from this feature to fit big models, i.e. 5000 individuals in 555 seconds = 9 minutes in a MacBook 4GB RAM. If you use this feature keep in mind that residuals and fitted values cannot be used given the roation of the data.}
    \item{forced}{a list of values for variance-covariance components to be used if the user wants to force those values. The format is the same than the init argument.}
   \item{complete}{a TRUE/FALSE statement to indicate if the function should impute the cases where at least for one trait there's an observation. Default is TRUE in order to use all data points where at least one observation exist across traits.}
   \item{restrained}{a numeric argument specifying which variance-covariance parameters should be restrained. Please see \code{\link{vctable}} and \code{\link{vctable.help}} functions if you want to use it although you shouldn't since is made by the \code{\link{mmer2}} function.}
   \item{na.method.X}{one of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the X matrices for fixed effects with the median value. If "exclude" is selected it will get rid of all rows with missing values for the X (fixed) covariates. The default is "exclude". The "include" option should be used carefully.}
   \item{na.method.Y}{one of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the Y matrix/dataframe with the median value. If "exclude" is selected it will get rid of missing values for the estimation of variance components. The default is "exclude".}
    \item{REML}{a TRUE/FALSE value to indicate if REML or ML should be used for optimization. Not functional yet. Only REML available.}
    
}
\details{

The package has been equiped with several datasets to learn how to use the sommer package: 

* \code{\link{HDdata}} and \code{\link{FDdata}} datasets have examples to fit half and full diallel designs. 

* \code{\link{h2}} to calculate heritability

* \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets to perform genomic prediction in hybrid single crosses 

* \code{\link{wheatLines}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link{CPdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. 
* \code{\link{PolyData}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link{gryphondata}} data contains an example of an animal model including pedigree information.

* \code{\link{BTdata}} dataset contains an animal (birds) model.

If you need to use pedigree you need to convert your pedigree into a relationship matrix (use the `getA` function from the pedigreemm package). Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the mmer-type of functions. Additional functions for genetic analysis have been included such as False Discovery Rate calculation (\code{\link{fdr}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}), phasing F1 or CP populations (\code{\link{phase.F1}}). We have also included the famous \code{\link{transp}} function to get transparent colors. Useful functions for analyzing such designs are included such as the \code{\link{blocker}} and \code{\link{fill.design}}. 

\strong{MODELS ENABLED}

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
\item{fish.inv}{inverse of the Fisher's information.}
\item{residuals}{Residual values e = Y - XB}
\item{cond.residuals}{Conditional residual e.cond = Y - (XB + ZU)}
\item{LL}{LogLikelihood}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{X}{ incidence matrix for fixed effects}
\item{dimos}{dimnensions for incidence matrix for random effects}
\item{sigma.scaled}{scaled variance-covariance components}
\item{fitted.y}{ Fitted values y.hat=XB+Zu}
\item{fitted.u}{ Fitted values only across random effects u.hat=Zu.1+....+Zu.i}
\item{ZETA}{Original incidence and variance covariance matrices used in the model fit.}
\item{K}{ variance-covariance matrix for random effects. If more than one random effect this will be the diagonal binding of individual K matrices.}
\item{method}{The method for extimation of variance components specified by the user.}
\item{forced}{a TRUE/FALSE statement indicating if user used the forced argument.}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{restrained}{table of restrained parameters.}


}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

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
#### breeding values for 3 variance component
#### with one variance covariance structure
####=========================================####
####=========================================####
data(CPdata)
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix 
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)
mix1 <- mmer2(Yield~1,
              random=~g(id) + Rowf + Colf, 
              rcov=~units,
              G=list(id=A), 
              data=CPpheno)
summary(mix1)
#### calculate heritability
pin(mix1, h1 ~ V1/(V1+V4) )
# assumes diag(trait) structure for univariate models

# #### for multivariate models
# mix2 <- mmer2(cbind(Yield,color)~1,
#               random = ~ us(trait):g(id) + 
#                             diag(trait):Rowf + 
#                                 diag(trait):Colf, 
#               rcov = ~ us(trait):units,
#               G=list(id=A), 
#               data=CPpheno)
# summary(mix2)

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

####=========================================####
####=========================================####
#### EXAMPLE 2
#### breeding values with 3 variance components
#### for hybrid prediction
####=========================================####
####=========================================####

# data(cornHybrid)
# hybrid2 <- cornHybrid$hybrid # extract cross data
# A <- cornHybrid$K
# K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)
# K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)
# S <- kronecker(K1, K2, make.dimnames=TRUE) ; dim(S)   
# 
# head(hybrid2)
# ans <- mmer2(Yield ~ Location, 
#              random = ~ g(GCA1) + g(GCA2) + g(SCA), 
#              rcov = ~ units,
#              G=list(GCA1=K1, GCA2=K2, SCA=S),
#              data=hybrid2)
# ans$var.comp
# summary(ans)
# 
# #### you can specify heterogeneus variances
# ans <- mmer2(Yield ~ 1, 
#              random = ~ GCA2 + at(Location):GCA2, 
#              rcov = ~ at(Location):units,
#              data=hybrid2)
# summary(ans)
# 
# ##### example of multivariate model
# ans <- mmer2(cbind(Yield,PlantHeight) ~ 1, 
#              random = ~ us(trait):GCA2 + us(trait):at(Location):GCA2, 
#              rcov = ~ diag(trait):at(Location):units,
#              data=hybrid2)
# summary(ans)

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

####=========================================####
#### EXAMPLE 3
#### different models with sommer
####=========================================####

# data(example)
# head(example)
# 
# #### Univariate homogeneous variance models ####
# 
# ## Compound simmetry (CS) model
# ans1 <- mmer2(Yield~Env, 
#               random= ~ Name + Env:Name,
#               rcov= ~ units,
#               data=example)
# summary(ans1)
# 
# #### Univariate heterogeneous variance models ####
# 
# ## Compound simmetry (CS) + Diagonal (DIAG) model
# ans2 <- mmer2(Yield~Env, 
#               random= ~Name + at(Env):Name,
#               rcov= ~ at(Env):units,
#               data=example)
# summary(ans2)
# 
# ##### Multivariate homogeneous variance models ####
# 
# ## Multivariate Compound simmetry (CS) model
# ans3 <- mmer2(cbind(Yield, Weight) ~ Env, 
#               random= ~ us(trait):Name + us(trait):Env:Name,
#               rcov= ~ us(trait):units,
#               data=example)
# summary(ans3)
# 
# ##### Multivariate heterogeneous variance models ####
# 
# ## Multivariate Compound simmetry (CS) + Diagonal (DIAG) model
# ans4 <- mmer2(cbind(Yield, Weight) ~ Env, 
#               random= ~ us(trait):Name + us(trait):at(Env):Name,
#               rcov= ~ us(trait):at(Env):units,
#               data=example)
# summary(ans4)

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

####=========================================####
#### EXAMPLE 4
#### comparison with lmer, install 'lme4' 
#### and run the code below
####=========================================####

#### lmer cannot use var-cov matrices so we will not 
#### use them in this comparison example

# library(lme4)
# library(sommer)
# data(cornHybrid)
# hybrid2 <- cornHybrid$hybrid
# fm1 <- lmer(Yield ~ Location + (1|GCA1) + (1|GCA2) + (1|SCA),
#             data=hybrid2 )
# out <- mmer2(Yield ~ Location, 
#              random = ~ GCA1 + GCA2 + SCA,
#              rcov = ~ units,
#              data=hybrid2)
# summary(fm1)
# summary(out)
# ### same BLUPs for GCA1, GCA2, SCA than lme4
# plot(out$cond.residuals, residuals(fm1))
# plot(out$u.hat[[1]], ranef(fm1)$GCA1[,1])
# plot(out$u.hat[[2]], ranef(fm1)$GCA2[,1])
# vv=which(abs(out$u.hat[[3]]) > 0)
# plot(out$u.hat[[3]][vv,], ranef(fm1)$SCA[,1])
# 
# ### a more complex model specifying which locations
# out2 <- mmer2(Yield ~ 1, 
#               random = ~ at(Location,c("3","4")):GCA2 + at(Location,c("3","4")):SCA,
#               rcov = ~ at(Location):units,
#               data=hybrid2)
# summary(out2)



}