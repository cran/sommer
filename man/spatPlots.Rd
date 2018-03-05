\name{spatPlots}
\alias{spatPlots}

\title{
Spatial plots
}
\description{
Plot fitted values for all terms in a model to assess spatial fit.
}
\usage{
spatPlots(object, by=NULL, colfunc=NULL,
         row="ROW",range="RANGE", wire=FALSE)

}

\arguments{
  \item{object}{a mixed model fitted using mmer2.}
\item{by}{
A character argument specifying the column name in the data frame in case the asreml model will specify a diagonal structure, meaning we will fit a different 2d spline at each FIELDINST
}
\item{colfunc}{
an optional function to creat colors for the spatial plots.
}
\item{row}{
name of the column in the dataset that identifies the spatial position in the x coordinate.
}
\item{range}{
name of the column in the dataset that identifies the spatial position in the y coordinate. 
}
\item{wire}{
a TRUE/FALSE statement indicating if the function should return a wire plot or the default levelplot.
}
}
\details{

The function only takes the asreml model and builds the fitted values for the terms required by building the incidence matrix and getting the BLUPs from the asreml model to build fitted = Zu. Then it makes a levelplot or wireplot to asses the spatial fit.

}
\value{
plots: spatial plots for all fitted values from the mmer2 model terms in the random part.

\describe{
\item{fits}{ A new dataset with the fitted values for all terms in the model.}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2017). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics (to appear). https://doi.org/10.1016/j.spasta.2017.10.003.

}

\examples{

data(CPdata)
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix 
A <- A.mat(CPgeno) # additive relationship matrix
#### look at the data and fit the model
head(CPpheno)

# mix1 <- mmer2(Yield~1,
#               random=~g(id)
#                       + Rowf + Colf
#                       + spl2D(Row,Col),
#               rcov=~units,
#               G=list(id=A),
#               data=CPpheno)
# summary(mix1)
# #### get the spatial plots
# spatPlots(mix1,row = "Row", range = "Col")

}
\seealso{\code{\link{mmer}} and \code{\link{mmer2}}-- the core functions of the package }
