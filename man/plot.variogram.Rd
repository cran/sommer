\name{plot.variogram.MMERM}
\alias{plot.variogram.MMERM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default variogram.MMERM plotting
}
\description{
Takes a fitted \code{variogram.MMERM} object produced by \code{variogram.MMERM()} and plots the associated sample variogram using an wireframe plot (package \code{lattice}).
}
\usage{
\method{plot}{variogram.MMERM}(x, stnd=TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{variogram.MMERM} as produced by \code{variogram.MMERM()}.}
  \item{stnd}{argument for ploting the residuals to know if they should be standarized.}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
This function as well as function \code{variogram.MMERM()} can only be used for regular two dimensional data.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Stefanova, K.T., Smith, A.B., and Cullis, B.R. (2009). Enhanced Diagnostics for the Spatial Analysis of Field Trials. Journal of Agricultural, Biological, and Environmental Statistics, 14, 392 - 410.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{mmer2}},\code{\link{mmer}}, \code{\link{variogram.MMERM}}
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
