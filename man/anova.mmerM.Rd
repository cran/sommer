\name{anova.mmerM}
\alias{anova.mmerM}
\title{anova form a GLMM fitted with mmer}
\description{\code{anova} method for class \code{"mmerM"}.
}
\usage{
\method{anova}{mmerM}(object, object2=NULL, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmerM"}}
 \item{object2}{an object of class \code{"mmerM"}, if NULL the program will return an error specifying that this function was developed for comparing mixed models using likelihood ratio tests and no to provide regular sum of squares results.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of anova
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\seealso{\code{\link{anova}}, \code{\link{mmer}}}
\keyword{models}

