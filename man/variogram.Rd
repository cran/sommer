\name{variogram}
\alias{variogram}
\alias{variogram.MMERM}
\title{
Sample variogram
}
\description{
Computes the sample variogram a dataset, usually residuals obtained from an mmer-type of functions.
}
\usage{

\method{variogram}{MMERM}(x, xcoor="ROW", ycoor="RANGE", 
     zcoor=NULL, by=NULL, ...)
}
\arguments{
  \item{x}{a data frame with the 3 dimensions (columns) to be plotted. Ideally for inspecting the residuals.}
  \item{xcoor}{a character argument indicating the name of the column in the dataset corresponding to the x coordinate.}
  \item{ycoor}{a character argument indicating the name of the column in the dataset corresponding to the y coordinate.}
  \item{zcoor}{a character argument indicating the name of the column in the dataset corresponding to the z coordinate.}
  \item{by}{a character argument indicating the name of the column in the dataset in case multiple residual plots are required and divided by a factor.}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented}
}
\details{
The present function computes the sample variogram on the basis of the (deviance) residuals of the fitted model. Currently, the function can only be applied for regular two-dimensional data, i.e, when the plots of the field are arranged in a regular two-dimensional array (usually defined by the column and row positions). 
For each pair of (deviance) residuals \eqn{e_i} and \eqn{e_j}, the half-squared difference is computed \deqn{v_{ij} = 0.5(e_i - e_j)^2,} as well as the corresponding column (\eqn{cd_{ij}}) and row displacements (\eqn{rd_{ij}}), with \deqn{cd_{ij} = |c_i - c_j|} and \deqn{rd_{ij} = |r_i - r_j|,} where \eqn{c_k} and \eqn{r_k} denote the column and row position of plot \eqn{k} respectively. The sample variogram is then defined as the triplet \deqn{(cd_{ij}, rd_{ij}, \bar{v}_{ij}),} where \eqn{\bar{v}_{ij}} denotes the average of the \eqn{v_{ij}} that share the same column and row displacements. 
For a more detailed description, see Gilmour et al. (1997).  
}
\value{
An object of class with the following components:
\item{data}{data frame including the following information: ``value'': the value of the sample variogram at each pair of column and row displacements; and ``length'': the number of observations used to compute the sample variogram at the corresponding pair of displacements.}
\item{col.displacement}{numerical vector containing the column displacements}
\item{row.displacement}{numerical vector containing the row displacements}
}
\references{
Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Stefanova, K.T., Smith, A.B. and Cullis, B.R. (2009). Enhanced Diagnostics for the Spatial Analysis of Field Trials. Journal of Agricultural, Biological, and Environmental Statistics, 14, 392 - 410.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

\code{\link{mmer2}}

}
\examples{

data(CPdata)
#### look at the data and fit the model
head(CPpheno)
mix1 <- mmer2(Yield~1,
              random=~ id
                      + Rowf + Colf,
              rcov=~units,
              data=CPpheno)
summary(mix1)

# compute the variogram
vm0 <- variogram(mix1, xcoor = "Row", ycoor = "Col")
# plot the variogram
plot(vm0$F1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
