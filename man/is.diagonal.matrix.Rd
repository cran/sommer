\name{is.diagonal.matrix}
\alias{is.diagonal.matrix}
\title{ Test for diagonal square matrix }
\description{
  This function returns \code{TRUE} if the given matrix argument x
  is a square numeric matrix and that the off-diagonal elements are close
  to zero in absolute value to within the given tolerance level. Otherwise,
  a \code{FALSE} value is returned.
}
\usage{
is.diagonal.matrix(x, tol = 1e-08)
}
\arguments{
  \item{x}{ a numeric square matrix }
  \item{tol}{ a numeric tolerance level usually left out }
}
\value{
  A TRUE or FALSE value.
}
\references{
  Bellman, R. (1987). \emph{Matrix Analysis}, Second edition, Classics in Applied Mathematics,
  Society for Industrial and Applied Mathematics.

  Horn, R. A. and C. R. Johnson (1990). \emph{Matrix Analysis}, Cambridge University Press.
}
\examples{
A <- diag( 1, 3 )
is.diagonal.matrix( A )
B <- matrix( c( 1, 2, 3, 4 ), nrow=2, byrow=TRUE )
is.diagonal.matrix( B )
C <- matrix( c( 1, 0, 0, 0 ), nrow=2, byrow=TRUE )
is.diagonal.matrix( C )
}
\keyword{ math }
