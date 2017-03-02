\name{yates.oats}
\alias{yates.oats}
\docType{data}
\title{
  Yield of oats in a split-block experiment
}
\description{
  The yield of oats from a split-plot field trial using three varieties
  and four levels of manurial treatment. The experiment was laid out in
  6 blocks of 3 main plots, each split into 4 sub-plots. The varieties
  were applied to the main plots and the manurial (nitrogen) treatments
  to the sub-plots.
}

\format{
  \describe{
    \item{\code{block}}{block factor with 6 levels}
    \item{\code{nitro}}{nitrogen treatment in hundredweight per acre}
    \item{\code{Variety}}{genotype factor, 3 levels}
    \item{\code{yield}}{yield in 1/4 lbs per sub-plot, each 1/80 acre.}
    \item{\code{row}}{row location}
    \item{\code{column}}{column location}
  }
}

\source{
  Yates, Frank (1935) Complex experiments,
  \emph{Journal of the Royal Statistical Society Suppl}. 2, 181--247.
}

\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S}. Fourth edition. Springer.
}

\examples{

### ========================== ###
### using the mmer2 function 
### ========================== ###

data(yates.oats)
head(yates.oats)
m3 <- mmer2(fixed=Y ~ V*N, random = ~ B + B:MP,
            data = yates.oats)
summary(m3)       
m3$var.comp

### ========================== ###
### using the mmer function 
### ========================== ###

###response
y<-yates.oats$Y
###fixed effects
X1 <- model.matrix(~V*N,yates.oats)
###random effects
Z1 <- model.matrix(~B-1, yates.oats)
Z2 <- model.matrix(~B:MP-1, yates.oats)
ETA <- list(B=list(Z=Z1),BMP=list(Z=Z2))
###run the model
m4 <- mmer(Y=y,X=X1,Z=ETA)
summary(m4)
}

