\name{hadamard.prod}
\alias{hadamard.prod}
\title{ Hadamard product of two matrices }
\description{
  This function returns the Hadamard or Shur product of two matrices, x and y, that have the same row and
  column dimensions.
}
\usage{
hadamard.prod(x, y)
}
\arguments{
  \item{x}{ a numeric matrix or vector object }
  \item{y}{ a numeric matrix or vector object }
}
\details{
  The Hadamard product is an element-by-element product of the two matrices. Let \eqn{{\bf{x}}} 
  and \eqn{{\bf{x}}} be two \eqn{m \times n} numeric matrices.  The Hadamard product is \eqn{{\bf{x}}\, \circ \,{\bf{y}} = \left[ {\begin{array}{*{20}{c}}
{{x_{1,1}}\,{y_{1,1}}}&{{x_{1,2}}\,{y_{1,2}}}& \cdots &{{x_{1,n}}\,{y_{1,n}}}\\
{{x_{2,1}}\,{y_{121}}}&{{x_{2,2}}\,{y_{2,2}}}& \cdots &{{x_{2,n}}\,{y_{2,n}}}\\
 \cdots & \cdots & \cdots & \cdots \\
{{x_{m,1}}\,{y_{m,1}}}&{{x_{m,2}}\,{y_{m,2}}}& \cdots &{{x_{m,n}}\,{y_{m,n}}}
\end{array}} \right]}.
  It uses the * operation in R.
}
\value{
  A matrix.
}
\references{
  Hadamard, J (1983). Resolution d'une question relative aux determinants, \emph{Bulletin des Sciences
  Mathematiques}, 17, 240-246.
  
  Styan, G. P. H. (1973). Hadamard Products and Multivariate Statistical Analysis,
  \emph{Linear Algebra and Its Applications}, Elsevier, 6, 217-240.
}
\note{
  The function converts vectors to matrices if necessary.
  The function stops running if x or y is not numeric and an error message is displayed.
  The function also stops running if x and y do not have the same row and column dimensions and an error mesage
  is displayed.
}
\examples{
x <- matrix( c( 1, 2, 3, 4 ), nrow=2, byrow=TRUE )
y <- matrix( c( 2, 4, 6, 8 ), nrow=2, byrow=TRUE )
z <- hadamard.prod( x, y )
print( z )
}
\keyword{ math }
