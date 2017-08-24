\name{residuals.MMERM}
\alias{residuals.MMERM}
\title{Residuals form a GLMM fitted with mmer}
\description{\code{residuals} method for class \code{"MMERM"}.
}
\usage{
\method{residuals}{MMERM}(object, type = "conditional", \dots)
}
\arguments{
 \item{object}{an object of class \code{"MMERM"}}
 \item{type}{ the type of residuals which should be returned. The alternatives are: "conditional" (e=y-(Xb+Zu)) and "regular" (e=y-Xb).}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\seealso{\code{\link{residuals}}, \code{\link{mmer}}}
\keyword{models}

