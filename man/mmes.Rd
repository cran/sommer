\name{mmes}
\alias{mmes}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations \strong{s}olver

}
\description{

The \code{mmes} function uses either the direct inversion or the Henderson mixed 
model equations algorithms coded in C++ using the Armadillo library to optimize matrix 
operations. For problems of the type  c > r (more coefficients to estimate than records available),
the direct inversion algorithm is faster (using the argument \code{henderson=FALSE}; default).
For more records than coefficients set \code{henderson=TRUE} and make sure that you provide the
relationship matrix as an inverse (see \code{\link{vsm}()} function for details).

}

\usage{

mmes(fixed, random, rcov, data, W, nIters=50, tolParConvLL = 1e-04,
     tolParConvNorm = 1e-04, tolParInv = 1e-06, naMethodX="exclude",
     naMethodY="exclude", returnParam=FALSE, dateWarning=TRUE,
     verbose=TRUE,addScaleParam=NULL, stepWeight=NULL, emWeight=NULL,
     contrasts=NULL, getPEV=TRUE, henderson=FALSE)
}
\arguments{
  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:

  \emph{response ~ covariate}

  }

  \item{random}{A formula specifying the name of the \strong{random effects}, e.g., 
  \emph{random= ~ genotype + year}.

  Useful functions can be used to fit heterogeneous variances and other special models
  (\emph{see 'Special Functions' in the Details section for more information}):

  \code{\link{vsm}(...,Gu)} is the main function to specify variance models and
  special structures for random effects. On the ... argument you provide the 
  unknown variance-covariance structures (e.g., usm,dsm,atm,csm) and the random 
  effect where such covariance structure will be used (the random effect of interest). 
  Gu is used to provide known covariance matrices among the levels of the random effect.
  Inverse matrices when the argument henderson=TRUE, and raw (non-inverse) matrices
  when henderson=FALSE (default direct inversion). Auxiliar functions for building
  the variance models are:

  ** \code{\link{dsm}(x)}, \code{\link{usm}(x)}, \code{\link{rrm}(x,y,z)} , 
  \code{\link{ism}(x)},\code{\link{csm}(x)}, and \code{\link{atm}(x,levs)} can 
  be used to specify unknown diagonal, unstructured, reduced-rank, identity, 
  and customized unstructured and diagonal covariance structures respectively 
  to be estimated by REML.

  ** \code{\link{unsm}(x)}, \code{\link{fixm}(x)} and \code{\link{diag}(x)} can
  be used to build easily matrices to specify constraints in the Gtc argument 
  of the \code{\link{vsm}()} function.

  ** \code{\link{overlay}()}, \code{\link{spl2Dc}()}, and \code{\link{leg}()},
  \code{\link{redmm}()} functions can be used to specify overlayed of design
  matrices of random effects, two dimensional spline, random regression, and 
  dimensionality-reduction models within the \code{\link{vsm}()} function.

  }
  \item{rcov}{A formula specifying the name of the \strong{error term}, e.g.,
  \emph{rcov= ~ units}.

  Special heterogeneous and special variance models and constraints for the 
  residual part are the same used on the random term but the name of the random 
  effect is always "units" which can be thought as a column with as many levels
  as rows in the data, e.g., \emph{rcov=~vsm(dsm(covariate),ism(units))}

  When fitting structures at the level of residuals please make sure that your 
  data is sorted based on the factors defining the structure. For example, for
  \emph{rcov= ~ vsm(dsm(xx), ism(units))} sort the datatset by the variable xx.

  }
  \item{data}{A data frame containing the variables specified in the formulas 
  for response, fixed, and random effects.}
  \item{W}{Weights matrix (e.g., when covariance among plots exist). Internally
  W is squared and inverted as Wsi = solve(chol(W)), then the residual matrix
  is calculated as R = Wsi*O*Wsi.t(), where * is the matrix product, and O 
  is the original residual matrix.}
   \item{nIters}{Maximum number of iterations allowed.}
  \item{tolParConvLL}{Convergence criteria based in the change of log-likelihood
  between iteration i and i-1.}
  \item{tolParConvNorm}{When using the Henderson method this argument is the 
  convergence criteria based in the norm proposed by Jensen, Madsen and
  Thompson (1997):

  e1 = || InfMatInv.diag()/sqrt(N) * dLu ||

  where InfMatInv.diag() is the diagonal of the inverse of the information
  matrix, N is the total number of variance components, and dLu is the vector 
  of first derivatives.

  }
  \item{tolParInv}{Tolerance parameter for matrix inverse used when singularities 
  are encountered in the estimation procedure. By default the value is 1e-06. This
  parameter should be fairly small because the it is used to bend matrices like the
  information matrix in the henderson algorithm or the coefficient matrix when it
  is not positive-definite.}
   \item{naMethodX}{One of the two possible values; "include" or "exclude". If 
   "include" is selected then the function will impute the X matrices for fixed
   effects with the median value. If "exclude" is selected it will get rid of 
   all rows with missing values for the X (fixed) covariates. The default is 
   "exclude". The "include" option should be used carefully.}
   \item{naMethodY}{One of the three possible values; "include", "include2" or 
   "exclude" (default) to treat the observations in response variable to be used 
   in the estimation of variance components. The first option "include" will 
   impute the response variables for all rows with the median value, whereas 
   "include2" imputes the responses only for rows where there is observation(s)
   for at least one of the responses (only available in the multi-response models). 
   If "exclude" is selected (default) it will get rid of rows in response(s) 
   where missing values are present for at least one of the responses.}
   \item{returnParam}{A TRUE/FALSE value to indicate if the program should 
   return the parameters to be used for fitting the model instead of fitting
   the model.}
   \item{dateWarning}{A TRUE/FALSE value to indicate if the program should warn 
   you when is time to update the sommer package.}
  \item{verbose}{A TRUE/FALSE value to indicate if the program should return 
  the progress of the iterative algorithm.}
  \item{addScaleParam}{additional scale parameters for the thetaF matrix when 
  using Henderson method (henderson=TRUE).}
  \item{stepWeight}{A vector of values (of length equal to the number of 
  iterations) indicating the weight used to multiply the update (delta) for
  variance components at each iteration. If NULL the 1st iteration will be 
  multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument can
  help to avoid that variance components go outside the parameter space in the 
  initial iterations which happens very often with the AI method but it can be 
  detected by looking at the behavior of the likelihood. In that case you may
  want to give a smaller weight.}
  \item{emWeight}{A vector of values (of length equal to the number of iterations) 
  indicating with values between 0 and 1 the weight assigned to the EM 
  information matrix. And the values 1 - emWeight will be applied to the AI
  information matrix to produce a joint information matrix. By default the 
  function gives a weight to the EM algorithm of a logarithmic decrease rate 
  using the following code:
  
  \code{ stan(logspace(seq(1,-1,- 2/nIters), p=3)) }. 
  
  }
  \item{contrasts}{an optional list. See the contrasts.arg of model.matrix.default.}
  \item{getPEV}{a logical value indicating if PEV should be returned when the 
  direct inversion algorithm is used. It does not apply when henderson argument is
  TRUE.}
  \item{henderson}{a logical value indicating if the solving algorithm should 
  be direct inversion (henderson is FALSE) or Henderson's method (henderson is
  TRUE). Default is direct inversion.}
}

\details{

The use of this function requires a good understanding of mixed models. 
Please review the 'sommer.quick.start' vignette and pay attention to details 
like format of your random and fixed variables (e.g. character and factor 
variables have different properties when returning BLUEs or BLUPs).

\strong{For tutorials} on how to perform different analysis with sommer please
look at the vignettes by typing in the terminal:

vignette("v1.sommer.quick.start")

vignette("v2.sommer.changes.and.faqs")

vignette("v3.sommer.qg")

vignette("v4.sommer.gxe")

\strong{Citation}

Type \emph{citation("sommer")} to know how to cite the sommer package in your publications.

\strong{Special variance structures}

 \code{\link{vsm}(\link{atm}(x,levels),ism(y))}

 can be used to specify heterogeneous variance for the "y" covariate at specific
 levels of the covariate "x", e.g., \emph{random=~vsm(at(Location,c("A","B")),ism(ID))}
 fits a variance component for ID at levels A and B of the covariate Location.

 \code{\link{vsm}(\link{dsm}(x),ism(y))}

 can be used to specify a diagonal covariance structure for the "y" covariate 
 for all levels of the covariate "x", e.g., \emph{random=~vsm(dsm(Location),ism(ID))} 
 fits a variance component for ID at all levels of the covariate Location.

 \code{\link{vsm}(\link{usm}(x),ism(y))}

 can be used to specify an unstructured covariance structure for the "y" covariate 
 for all levels of the  covariate "x", e.g., \emph{random=~vsm(usm(Location),ism(ID))}  
 fits variance and covariance components for ID at all levels of the covariate Location.
 
 \code{\link{vsm}(\link{usm}(\link{rrm}(x,y,z,nPC)),ism(y))}

 can be used to specify an unstructured covariance structure for the "y" effect
 for all levels of the  covariate "x", and a response variable "z", e.g., 
 \emph{random=~vsm(rrm(Location,ID,response, nPC=2),ism(ID))}  fits a reduced-rank 
 factor analytic covariance for ID at 2 principal components of the covariate Location.

  \code{\link{vsm}(ism(\link{overlay}(...,rlist=NULL,prefix=NULL)))}

 can be used to specify overlay of design matrices between consecutive random 
 effects specified, e.g., \emph{random=~vsm(ism(overlay(male,female)))} overlays
 (overlaps) the incidence matrices for the male and female random effects to obtain
 a single variance component for both effects. The `rlist` argument is a list with 
 each element being a numeric value that multiplies the incidence matrix to be 
 overlayed. See \code{\link{overlay}} for details.Can be combined with vsm().
 
   \code{\link{vsm}(ism(\link{redmm}(x,M,nPC)))}

 can be used to create a reduced model matrix of an effect (x) assumed to be a 
 linear function of some feature matrix (M), e.g., \emph{random=~vsm(ism(redmm(x,M)))} 
 creates an incidence matrix from a very large set of features (M) that belong 
 to the levels of x to create a reduced model matrix. See \code{\link{redmm}} 
 for details.Can be combined with vsm().

  \code{\link{vsm}(\link{leg}(x,n),ism(y))}

  can be used to fit a random regression model using a numerical variable \code{x} 
  that marks the trayectory for the random effect \code{y}. The leg function can 
  be combined with the special functions \code{dsm}, \code{usm} \code{at} and \code{csm}. 
  For example \emph{random=~vsm(leg(x,1),ism(y))} or \emph{random=~vsm(usm(leg(x,1)),ism(y))}.

  \code{\link{spl2Dc}(x.coord, y.coord, at.var, at.levels))}

  can be used to fit a 2-dimensional spline (e.g., spatial modeling) using 
  coordinates \code{x.coord} and \code{y.coord} (in numeric class) assuming 
  multiple variance components. The 2D spline can be fitted at specific 
  levels using the \code{at.var} and \code{at.levels} arguments. For example
  \emph{random=~spl2Dc(x.coord=Row,y.coord=Range,at.var=FIELD)}.

\strong{Covariance between random effects}

 \code{\link{covm}( \link{vsm}(\link{ism}(ran1)), \link{vsm}(\link{ism}(ran2)) )}

 can be used to specify covariance between two different random effects, e.g.,
 \emph{random=~covm( vsm(ism(x1)), vsm(ism(x2)) )} where two random effects in 
 their own vsm() structure are encapsulated. Only applies for simple random effects.
 
\strong{S3 methods}

S3 methods are available for some parameter extraction such as
\code{\link{fitted.mmes}}, \code{\link{residuals.mmes}}, \code{\link{summary.mmes}},
\code{\link{randef}}, \code{\link{coef.mmes}}, \code{\link{anova.mmes}},
\code{\link{plot.mmes}}, and \code{\link{predict.mmes}} to obtain adjusted means.
In addition, the \code{\link{vpredict}} function (replacement of the pin 
function) can be used to estimate standard errors for linear combinations of 
variance components (e.g., ratios like h2). The \code{\link{r2}} function calculates
reliability. 


\strong{Additional Functions}

Additional functions for genetic analysis have been included such as relationship 
matrix building (\code{\link{A.mat}}, \code{\link{D.mat}}, \code{\link{E.mat}},
\code{\link{H.mat}}),  build a genotypic hybrid marker matrix (\code{\link{build.HMM}}), 
plot of genetic maps (\code{\link{map.plot}}), and manhattan plots
(\code{\link{manhattan}}). If you need to build a pedigree-based relationship 
matrix use the \code{getA} function from the pedigreemm package.

\strong{Bug report and contact}

If you have any technical questions or suggestions please post it in
https://stackoverflow.com or https://stats.stackexchange.com

If you have any bug report please go to https://github.com/covaruber/sommer or 
send me an email to address it asap, just make sure you have read the vignettes 
carefully before sending your question.

\strong{Example Datasets}

The package has been equiped with several datasets to learn how to use the sommer package:

* \code{\link{DT_halfdiallel}}, \code{\link{DT_fulldiallel}} and \code{\link{DT_mohring}} 
datasets have examples to fit half and full diallel designs.

* \code{\link{DT_h2}} to calculate heritability

* \code{\link{DT_cornhybrids}} and \code{\link{DT_technow}} datasets to perform
genomic prediction in hybrid single crosses

* \code{\link{DT_wheat}} dataset to do genomic prediction in single crosses in 
species displaying only additive effects.

* \code{\link{DT_cpdata}} dataset to fit genomic prediction models within a
biparental population coming from 2 highly heterozygous parents including additive, 
dominance and epistatic effects.

* \code{\link{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids.

* \code{\link{DT_gryphon}} data contains an example of an animal model including 
pedigree information.

* \code{\link{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link{DT_legendre}} simulated dataset for random regression model.

* \code{\link{DT_sleepstudy}} dataset to know how to translate lme4 models to 
sommer models.

* \code{\link{DT_ige}} dataset to show how to fit indirect genetic effect models.

\strong{Models Enabled}

For details about the models enabled and more information about the covariance 
structures please check the help page of the package (\code{\link{sommer}}).

}

\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{data}{the dataset used in the model fitting.}
\item{Dtable}{the table to be used for the predict function to help the program
              recognize the factors available.}
\item{llik}{the vector of log-likelihoods across iterations}
\item{b}{the vector of fixed effect.}
\item{u}{the vector of random effect.}
\item{bu}{the vector of fixed and random effects together.}
\item{Ci}{the inverse of the coefficient matrix.}
\item{Ci_11}{the inverse of the coefficient matrix pertaining to the fixed effects.}
\item{theta}{a list of estimated variance covariance matrices. Each element of 
             the list corresponds to the different random and residual components}
\item{theta_se}{inverse of the information matrix.}
\item{InfMat}{information matrix.}
\item{monitor}{The values of the variance-covariance components across iterations 
              during the REML estimation.}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{partitions}{a list where each element contains a matrix indicating where 
                  each random effect starts and ends.}
\item{partitionsX}{a list where each element contains a matrix indicating where 
                   each fixed effect starts and ends.}
\item{percDelta}{the matrix of percentage change in deltas (see tolParConvNorm argument).}
\item{normMonitor}{the matrix of the three norms calculated (see tolParConvNorm argument).}
\item{toBoundary}{the matrix of variance components that were forced to the
                   boundary across iterations.}
\item{Cchol}{the Cholesky decomposition of the coefficient matrix.}
\item{y}{the response vector.}
\item{W}{the column binded matrix W = [X Z] }
\item{uList}{a list containing the BLUPs in data frame format where rows are 
             levels of the random effects and column the different factors at 
             which the random effect is fitted. This is specially useful for 
             diagonal and unstructured models.}
\item{uPevList}{a list containing the BLUPs in data frame format where rows are 
                levels of the random effects and column the different factors at 
                which the random effect is fitted. This is specially useful for 
                diagonal and unstructured models.}
\item{args}{the fixed, random and residual formulas from the mmes model.}
\item{constraints}{The vector of constraints.}

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using 
the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744

Jensen, J., Mantysaari, E. A., Madsen, P., and Thompson, R. (1997). Residual 
maximum likelihood estimation of (co) variance components in multivariate mixed 
linear models using average information. Journal of the Indian Society of
Agricultural Statistics, 49, 215-236.

Sanderson, C., & Curtin, R. (2025). Armadillo: An Efficient Framework for 
Numerical Linear Algebra. arXiv preprint arXiv:2502.03000.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for 
variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

}
\author{

Coded by Giovanny Covarrubias-Pazaran with contributions of Christelle Fernandez 
Camacho to the henderson algorithm.

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####

data(DT_example)
DT <- DT_example
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)

# ####===========================================####
# #### Univariate heterogeneous variance models  ####
# ####===========================================####
# DT=DT[with(DT, order(Env)), ]
# ## Compound simmetry (CS) + Diagonal (DIAG) model
# ans2 <- mmes(Yield~Env,
#              random= ~Name + vsm(dsm(Env),ism(Name)),
#              rcov= ~ vsm(dsm(Env),ism(units)),
#              data=DT)
# summary(ans2)
# 
# ####===========================================####
# ####  Univariate unstructured variance models  ####
# ####===========================================####
# 
# ans3 <- mmes(Yield~Env,
#              random=~ vsm(usm(Env),ism(Name)),
#              rcov=~vsm(dsm(Env),ism(units)),
#              data=DT)
# summary(ans3)

}
