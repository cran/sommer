\name{mmer2}
\alias{mmer2}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations in \strong{R}

}
\description{

Fits a multivariate or univariate linear mixed model by likelihood methods (REML). It has been \strong{implemented to work with a data frame in a formula-based fashion} which will create the incidence and variance covariance matrices for the user and call the \code{\link{mmer}} function in the background. The optimization methods are; Direct-Inversion Newton-Raphson (NR) (Tunnicliffe 1989; Lee et al. 2015; Maier et al. 2015), and Efficient Mixed Model Association (EMMA) (Kang et al. 2008). These algorithms are \strong{intended to be used for problems of the type p > n}. For problems with sparse covariance structures, or problems of the type n > p, the MME-based algorithms are faster and we recommend to shift to the use of MME-based software (i.e. lme4, breedR, or asreml-R).

The package also provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices to model covariances among genotypes for plant and animal breeding problems.Spatial modeling can be made through the two dimensional splines (see details below).

\strong{For tutorials} on how to perform different analysis with sommer please look at the vignettes by typing in the terminal:

\strong{vignette("sommer.start")} 

\strong{vignette("sommer")}

}
\usage{

mmer2(fixed, random, rcov, data, weights, G=NULL, 
      grouping=NULL, method="NR", init=NULL, iters=20,
      tolpar = 1e-06, tolparinv = 1e-06, draw=FALSE, 
      silent=FALSE, constraint=TRUE, forced=NULL, 
      complete=TRUE, restrained=NULL,na.method.X="exclude",
      na.method.Y="exclude", REML=TRUE, init.equal=TRUE,
      return.param=FALSE,date.warning=TRUE)
}
\arguments{
  \item{fixed}{a formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:
  
  \emph{Yield ~ Location} for univariate models
  
  \emph{cbind(Yield,color) ~ Location} for multivariate models.}
  
  \item{random}{a formula specifying the name of the \strong{random effects}, i.e. \emph{random= ~ genotype + year}.
  
  Useful functions can be used to fit heterogeneous variances and other special models (\emph{see 'Special Functions' in the Details section for more information}):
  
  \code{diag(trait)} and \code{us(trait)} can be used to specify diagonal and unstructured covariance structures among responses in multivariate models.
  
  \code{at}(.) can be used to specify heterogeneous variance for specific levels of a random effect
  
  \code{diag}(.), \code{us}(.) can be used to specify diagonal and unstructured covariance structures among levels of a random effect.
  
  \code{and}(.)   can be used to specify overlay of design matrices of random effects
  
  \code{g}(.)   can be used to specify that a known covariance structure for a random effect
  
  \code{grp}(.)   can be used to specify customized incidence matrices for random effects
  
  \code{spl2D}(.)   can be used to fit a 2-dimensional spline (i.e. spatial modeling; see Special functions section below).
  
  }
  \item{rcov}{a formula specifying the name of the \strong{error term}, i.e. \emph{rcov= ~ units} or \emph{rcov=~at(.):units}.
  
  The functions that can be used to fit heterogeneous residual variances are:
  
    \code{diag(trait)} and \code{us(trait)} can be used to specify diagonal and unstructured covariance structures among responses for the residuals in multivariate models.
  
  \code{at}(.)   can be used to specify heterogeneous variance, i.e. \emph{rcov=~at(Location):units}
  
  }
  \item{data}{a data frame containing the variables specified in the formulas for response, fixed, and random effects.}
  \item{weights}{name of the covariate for weights. To be used for the product R = Wsi*R*Wsi, where * is the matrix product, Wsi is the square root of the inverse of W and R is the residual matrix.}
  \item{G}{a list containing the \strong{variance-covariance matrices for the random effects (not inverses)}, i.e. if the user specifies: 
  
  \emph{random=~ g(genotype) + g(year)}
  
  then the G argument should be used as:  
  
  \emph{G=list(genotype=A1, year=A2)} 
  
  where A1 and A2 are the variance-covariance matrices for random effects 'genotype' and 'year' respectively. Opposite to asreml you provide the original var-covar matrix, not its inverse.
  
  For special cases such as the \code{overlay()} and \code{grp()} functions. The variance covariance matrices in G should be provided with the entire name and the random formula should NOT use the \code{g()} function for such random effect, i.e:
  
  \emph{random=~overlay(females, males), G=list('overlay(females, males)'=A3)}
  
  }
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
    \item{forced}{a list of values for variance-covariance components to be used if the user wants to force those values. The format is the same than the init argument.}
   \item{complete}{a TRUE/FALSE statement to indicate if the function should impute the cases where at least for one trait there's an observation. Default is TRUE in order to use all data points where at least one observation exist across traits.}
   \item{restrained}{a numeric argument specifying which variance-covariance parameters should be restrained. Please see \code{\link{vctable}} and \code{\link{vctable.help}} functions if you want to use it although you shouldn't since is made by the \code{\link{mmer2}} function.}
   \item{na.method.X}{one of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the X matrices for fixed effects with the median value. If "exclude" is selected it will get rid of all rows with missing values for the X (fixed) covariates. The default is "exclude". The "include" option should be used carefully.}
   \item{na.method.Y}{one of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the Y matrix/dataframe with the median value. If "exclude" is selected it will get rid of missing values for the estimation of variance components. The default is "exclude".}
    \item{REML}{a TRUE/FALSE value to indicate if REML or ML should be used for optimization. Not functional yet. Only REML available.}
   \item{init.equal}{a TRUE/FALSE value to indicate if the program should use the same initial values for all variance-covariance components in the multivariate models. Default is TRUE, otherwise the function will calculate the variance-covariance with the raw data and use them as initial values.}
   \item{return.param}{a TRUE/FALSE value to indicate if the program should return the parameters used for modeling without fitting the model.}
   \item{date.warning}{a TRUE/FALSE value to indicate if the program should warn you when is time to update the sommer package.}
    
}

\details{

\strong{Special Functions}

 \code{at(x,levels):y}   
 
 can be used to specify heterogeneous variance for the "y"" covariate at specific levels of the factor covariate "x", i.e. \emph{random=~at(Location,c("A","B")):ID} fits a variance component for ID at levels A and B of the factor covariate Location.
  
 \code{diag(x):y}   
 
 can be used to specify a diagonal covariance structure for the "y"" covariate for all levels of the factor covariate "x", i.e. \emph{random=~diag(Location):ID}  fits a variance component for ID at all levels of the factor covariate Location.
 
  \code{diag(trait):x}  \strong{or}   \code{diag(trait):f(x):y}
 
 can be used to specify a diagonal covariance structure among responses in multivariate models for the random effect "x" or "f(x):y" where f is an additional function applied at the random effect level (i.e. diag(x):y), i.e. \emph{random=~diag(trait):ID} along with a multivariate response fits a diagonal covariance structure among responses for the random effect ID.
  
 \code{us(x):y}   
 
 can be used to specify an unstructured covariance structure for the "y"" covariate for all levels of the factor covariate "x", i.e. \emph{random=~us(Location):ID}  fits variance and covariance components for ID at all levels of the factor covariate Location.
 
  \code{us(trait):x} \strong{or} \code{us(trait):f(x):y}
 
 can be used to specify an unstructured covariance structure among responses in multivariate models for the random effect "x" or "f(x):y" where f is an additional function applied at the random effect level (i.e. us(x):y),, i.e. \emph{random=~us(trait):ID} along with a multivariate response fits an unstructured covariance structure among responses for the random effect ID.
 
  \code{overlay(...,rlist=NULL,prefix=NULL)}   
 
 can be used to specify overlay of design matrices between consecutive random effects specified, i.e. \emph{random=~overlay(male,female)} overlays (overlaps) the incidence matrices for the male and female random effects to obtain a single variance component for both effects. The `rlist` argument is a list with each element being a numeric value that multiplies the incidence matrix to be overlayed. See \code{\link{overlay}} for details.
  
  \code{g(x)}   
 
 can be used to specify that the random effect "x" has a known covariance structure that will be provided in the \code{G} argument, i.e. \emph{random=~g(male), G=list(male=Am)} would indicate that the levels of the random effect called "male"" have a known covariance structure and is provided as a matrix called "Am" in the \code{G} argument.
 
  \code{eig(g(x))}   
 
 can be used to specify that we should apply the eigen decomposition to known covariance structure of the random effect "x" that was provided in the \code{G} argument, i.e. \emph{random=~ eig(g(id)), G=list(id=Am)} would perform the eigen decomposition proposed by Lee et al. (2015).
 
  \code{grp(x)}   
 
 can be used to specify customized incidence matrices for a random effect called "x" that will be provided in the \code{grouping} argument, i.e.  i.e. \emph{random=~ grp(ID), grouping=list(re=Zre)} would indicate that a new random effect called ID will be provided by the user an incidence matrix (i.e. Zre) in the \code{grouping} argument.
  
  \code{spl2D(x.coord, y.coord, at, at.levels, type, 
      nseg, pord, degree, 
      nest.div)}
      
  can be used to fit a 2-dimensional spline (i.e. spatial modeling) using coordinates \code{x.coord} and \code{y.coord}. The 2D spline can be fitted at specific levels using the \code{at} and \code{at.levels} arguments. For example \emph{random=~spl2D(x.coord=Row,y.coord=Range,at=FIELD)}.Two types of P-Spline models are available; "SAP" or "PSANOVA" using the \code{type} argument. The \code{nseg} argument controls the number of segments for each marginal (strictly nseg - 1 is the number of internal knots in the domain of the covariate), and \code{pord} the penalty order for each marginal. \code{degree} is the order of the polynomial of the B-spline basis for each marginal being the default 3 (cubic B-splines). The \code{nest.div} argument refers to the number of segments to be used for the construction of the nested B-spline basis for the smooth-by-smooth interaction component. Default set to 1, which implies that nested basis are not used. \code{nseg}, \code{pord}, \code{degree}, \code{nest.div} are vectors of length 2 to be applied to \code{x.coord} and \code{y.coord} respectively. See an example in the \code{\link{spatPlots}} documentation.

For a short tutorial on how to use this special functions you can look at the vignettes by typing in the terminal:

\emph{vignette('sommer.start')}

\strong{Example Datasets}

The package has been equiped with several datasets to learn how to use the sommer package: 

* \code{\link{HDdata}} and \code{\link{FDdata}} datasets have examples to fit half and full diallel designs. 

* \code{\link{h2example}} to calculate heritability

* \code{\link{cornHybrid}} and \code{\link{Technow_data}} datasets to perform genomic prediction in hybrid single crosses 

* \code{\link{wheatLines}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link{CPdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects.
* \code{\link{PolyData}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link{gryphondata}} data contains an example of an animal model including pedigree information.

* \code{\link{BTdata}} dataset contains an animal (birds) information.

* \code{\link{ExpDesigns}} dataset contains examples for fitting different experimental designs.

* \code{\link{spatPlots}} dataset contains examples of spatial modeling of field trials through the two-dimensional splines.

\strong{Additional Functions}

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{variogram}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the mmer-type of functions. 

Additional functions for genetic analysis have been included such as heritability (\code{\link{h2.fun}}), False Discovery Rate calculation (\code{\link{fdr}}), build a genotypic hybrid marker matrix (\code{\link{build.HMM}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}). If you need to use pedigree you need to convert your pedigree into a relationship matrix (i.e. use the `getA` function from the pedigreemm package). 

Useful functions for analyzing field trials are included such as the \code{spl2D}, \code{\link{spatPlots}}, \code{\link{blockerT}}, \code{\link{blockerL}}, and \code{\link{fill.design}}.

\strong{Models Enabled}

For details about the models enabled and more information about the covariance structures please check the help page ot the package (\code{\link{sommer}}).

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
\item{fish.inv}{inverse of the Fisher's information scaled.}
\item{fish.inv.nonscale}{inverse of the Fisher's information without scaling.}
\item{residuals}{Residual values e = Y - XB}
\item{cond.residuals}{Conditional residual e.cond = Y - (XB + ZU)}
\item{LL}{LogLikelihood}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{Y}{ response matrix}
\item{X}{ incidence matrix for fixed effects}
\item{dimos}{dimnensions for incidence matrix for random effects}
\item{sigma.scaled}{scaled variance-covariance components}
\item{sigma}{vectorized variance-covariance components}
\item{fitted.y}{ Fitted values y.hat=XB+Zu}
\item{fitted.u}{ Fitted values only across random effects u.hat=Zu.1+....+Zu.i}
\item{Zus}{Fitted values for each random effect by separate.}
\item{ZETA}{Original incidence and variance covariance matrices used in the model fit.}
\item{method}{The method for extimation of variance components specified by the user.}
\item{forced}{a TRUE/FALSE statement indicating if user used the forced argument.}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{restrained}{table of restrained parameters.}
\item{call}{formula for fixed, random and rcov used.}
\item{data}{original dataset provided.}


}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

Rodriguez-Alvarez, Maria Xose, et al. Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics 23 (2018): 52-71.

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