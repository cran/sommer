\name{poe}
\alias{poe}
\title{Short poems from Latin America, and other places why not?}
\description{
This function was designed for fun to remember the scientist (and me) that life is more than analyzing data. If you want to contribute with a poem you like for the following versions of sommer sent it to me to:

covarrubiasp@wisc.edu

}
\usage{
poe(h)
}

\arguments{
  \item{h}{a numeric value between 1 to 8 to select the poem. In the function a random sample is taken}
}
\value{
\describe{

\item{$poe.value}{ a poem to remember the scientist (and me) that life is more than analyzing data.}
}
}
\references{

Many

}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
poe(sample(1:8,1))
}