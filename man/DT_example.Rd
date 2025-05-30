\name{DT_example}
\alias{DT_example}
\alias{A_example}
\docType{data}
\title{
Broad sense heritability calculation.
}
\description{
This dataset contains phenotpic data for 41 potato lines evaluated in 3 environments in an RCBD design. The phenotypic trait is tuber quality and we show how to obtain an estimate of DT_example for the trait.

}
\usage{data("DT_example")}
\format{
  The format is:
 chr "DT_example"
}
\source{
This data was generated by a potato study.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
####=========================================####
#### EXAMPLES
#### Different models with sommer
####=========================================####

data(DT_example)
DT <- DT_example
A <- A_example
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)$varcomp

# ####===========================================####
# #### Univariate heterogeneous variance models  ####
# ####===========================================####
# 
# ## Compound simmetry (CS) + Diagonal (DIAG) model
# DT=DT[with(DT, order(Env)), ]
# ans2 <- mmes(Yield~Env,
#              random= ~Name + vsm(dsm(Env),ism(Name)),
#              rcov= ~ vsm(dsm(Env),ism(units)),
#              data=DT)
# summary(ans2)

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmes}}  }
