\name{is.square.matrix}
\alias{is.square.matrix}
\title{ Test for square matrix }
\description{
  The function returns TRUE if the argument is a square matrix and FALSE otherwise.
}
\usage{
is.square.matrix(x)
}
\arguments{
  \item{x}{ a matrix }
}
\value{
  TRUE or FALSE
}
\references{
  Bellman, R. (1987). \emph{Matrix Analysis}, Second edition, Classics in Applied Mathematics,
  Society for Industrial and Applied Mathematics.
}
\examples{
A <- matrix( seq( 1, 12, 1 ), nrow=3, byrow=TRUE )
is.square.matrix( A )
B <- matrix( seq( 1, 16, 1 ), nrow=4, byrow=TRUE )
is.square.matrix( B )
}
\keyword{ math }
