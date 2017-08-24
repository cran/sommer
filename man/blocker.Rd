\name{blocker}
\alias{blocker}
\title{Applying postblocking in a field}
\description{
  The blocker function takes a dataframe that contains columns indicating x and y coordinates and use them to create transversal blocking with different numbers or row width. For example, a square field with no blocks after applying the function will have 2 new columns (forward and backwards) for the transversal blocking desired (angle degree and #rows specified).
}
\usage{
  blocker(dat,rows="ROW",ranges="RANGE",by=NULL)
}

\arguments{
    \item{dat}{a dataframe with 2 obligatory columns; rows and ranges which can have different names and can be matched with the following arguments.}
  \item{rows}{the name of the numeric column that indicates the x coordinate direction in the field.}
  \item{ranges}{the name of the numeric column that indicates the y coordinate other direction in the field.}
  \item{by}{optional argument to indicate the name of the column of the dataframe x that indicates the environments so the field is filled by environment.}
  
}
\value{
  \describe{
    
    \item{$fin}{ a new dataframe with new columns for 2-direction blocking.}
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
    Fikret Isik. 2009. Analysis of Diallel Mating Designs. North Carolina State University, Raleigh, USA.
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(CPdata)
#### look at the data
head(CPpheno)
#### fill the design
gg <- fill.design(x=CPpheno, rows="Row",ranges="Col")
head(gg)
#### apply the postblocking 
gg2 <- blocker(dat=gg, rows="Col",ranges="Row")
head(gg2)
#### see the new blocking
# lattice::levelplot(block.1~Row*Col|FIELDINST, data=gg2)
# lattice::levelplot(block.2~Row*Col|FIELDINST, data=gg2)
# lattice::levelplot(block.3~Row*Col|FIELDINST, data=gg2)
#### now you can use them in your mixed models


}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}