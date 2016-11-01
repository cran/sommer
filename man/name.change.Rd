\name{name.change}
\alias{name.change}
\title{renaming a vector by adding zeros}
\description{
This function is a very simple function to rename a vector that was named wrong and zeros were not added by mistake, i.e.:

x <- c("P1","P9","P10","P99","P187")

but the user desires the a vector of the form:

x <- c("P001","P009","P010","P099","P187")

this function does that simple task.

}
\usage{
name.change(x, separ="P", maxn=999)
}

\arguments{
  \item{x}{a character vector with the names.}
  \item{separ}{a character argument indicating what comes before the numbers that need to be corrected.}
  \item{maxn}{maximum number of individual-label found in the vector. This will be important for the function to know how many zeros should be added.}
}
\value{
\describe{

\item{$x}{a character vector with the names corrected by the added zeros.}

}
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
x <- c("P1","P9","P10","P99","P187")
## try adding zeros by indicating that 999 could be the max number
name.change(x, sep="P",maxn=999)
name.change(x, sep="P",maxn=9999)
}
\seealso{The core function of the package \code{\link{mmer}}}