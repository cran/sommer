\name{fitted.mmec}
\alias{fitted.mmec}
\title{fitted form a LMM fitted with mmec}
\description{\code{fitted} method for class \code{"mmec"}.
}
\usage{
\method{fitted}{mmec}(object, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmec"}}
 \item{\dots}{Further arguments to be passed to the mmec function}
}
\value{
  vector of fitted values of the form y.hat = Xb + Zu including all terms of the model.
 }
\author{Giovanny Covarrubias}
\examples{
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmer(Yield~1,
#               random=~vsr(id,Gu=A)
#                       + Rowf + Colf + spl2Da(Row,Col),
#               rcov=~units,
#               data=DT)
# 
# ff=fitted(mix1)
# 
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# lattice::wireframe(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted,  
#           aspect=c(61/87,0.4), drape=TRUE,# col.regions = colfunc,
#           light.source=c(10,0,10))
# lattice::levelplot(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted, col.regions = colfunc)

}
\seealso{\code{\link{fitted}}, \code{\link{mmec}}}
\keyword{models}

