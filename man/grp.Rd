\name{grp}
\alias{grp}
\title{grouping functionality}
\description{
  This function allows mmer2 to recognize that an extra random effect not present in the dataset but provided by the user should be included in the model. The user specifies the name of the random effect within the grp function and provides the incidence matrix in the 'grouping' argument and its covariance structure if other than diagonal in the 'G' argument as any other random effect.
}
\usage{
  grp(x)
}

\arguments{
  \item{x}{Name of the random effect that will be provided in the grouping argument.}
}
\value{
  \describe{
    
    \item{$dd}{ an incidence matrix with columns for the specific levels specified.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{
data(CPdata)
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix 
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)
#### create a customized incidence matrix
Zr <- model.matrix(~Rowf-1, CPpheno); dim(Zr)
Kr <- diag(ncol(Zr));
colnames(Kr) <- rownames(Kr) <- colnames(Zr)

mix1 <- mmer2(Yield~1,
              random=~ Colf + grp(Rowf), 
              rcov=~units,
              grouping=list(Rowf=Zr),
              G=list(Rowf=Kr),
              data=CPpheno)
summary(mix1)

# mix2 <- mmer2(cbind(color,Yield)~1,
#               random=~diag(trait):Colf + us(trait):grp(g(Rowf)), 
#               rcov=~us(trait):units,
#               grouping=list(Rowf=Zr),
#               G=list(Rowf=Kr),
#               data=CPpheno)
# summary(mix2)
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}} }