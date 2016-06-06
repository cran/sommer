\name{matrix.trace}
\alias{matrix.trace}
\title{ The trace of a matrix }
\description{
  This function returns the trace of a given square numeric matrix.
}
\usage{
matrix.trace(x)
}
\arguments{
  \item{x}{ a matrix }
}
\value{
  A numeric value which is the sum of the values on the diagonal.
}
\references{
  Bellman, R. (1987). \emph{Matrix Analysis}, Second edition, Classics in Applied Mathematics,
  Society for Industrial and Applied Mathematics.
}
\note{
  If the argument x is not numeric, the function presents and error message and terminates.
  If the argument x is not a square matrix, the function presents an error message and terminates.
}
\examples{
A <- matrix( seq( 1, 16, 1 ), nrow=4, byrow=TRUE )
matrix.trace( A )
}
\keyword{ math }
