\name{residuals.mmerM}
\alias{residuals.mmerM}
\title{Residuals form a GLMM fitted with mmer}
\description{\code{residuals} method for class \code{"mmerM"}.
}
\usage{
\method{residuals}{mmerM}(object, type = "conditional", \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmerM"}}
 \item{type}{ the type of residuals which should be returned. The alternatives are: "conditional" (e=y-(Xb+Zu)) and "regular" (e=y-Xb).}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\seealso{\code{\link{residuals}}, \code{\link{mmer}}}
\keyword{models}

