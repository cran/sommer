\name{brewer.pal}
\alias{brewer.pal}
\title{Generate a sequence of colors for groups.}
\usage{
  brewer.pal(n, name)
}
\arguments{
  \item{n}{Number of different colors in the palette, minimum 3, maximum depending on palette.}

  \item{name}{A palette name from the lists below:
  Accent 8
Dark2 8
Paired 12
Pastel1 9
Pastel2 8
Set1 9
Set2 8
Set3 12
  }
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

