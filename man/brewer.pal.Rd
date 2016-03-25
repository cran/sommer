\name{brewer.pal}
\alias{brewer.pal}
\title{Generate a sequence of colors for groups.}
\usage{
  brewer.pal(n, name)
}
\arguments{
  \item{n}{Number of different colors in the palette, minimum 3, maximum depending on palette.}

  \item{name}{A palette name from the lists below.}
}
\value{
  A vector of colors.
}
\description{
  wrapper of brewer.pal function from RColorBrewer.
}
\examples{
{
# Plot a colorbar with brewer.pal
mypalette<-brewer.pal(7,"Greens")
}
}

