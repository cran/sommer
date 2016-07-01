\name{fitted.MMERM}
\alias{fitted.MMERM}
\title{fitted form a GLMM fitted with mmer}
\description{\code{fitted} method for class \code{"MMERM"}.
}
\usage{
\method{fitted}{MMERM}(object, type = "complete", \dots)
}
\arguments{
 \item{object}{an object of class \code{"MMERM"}}
 \item{type}{ the type of fitted which should be returned. The alternatives are: "complete" (y.hat=Xb+Zu) and "incomplete" (y.hat=Zu).}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of fitted
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\seealso{\code{\link{fitted}}, \code{\link{mmer}}}
\keyword{models}

