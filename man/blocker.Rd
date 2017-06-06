\name{blocker}
\alias{blocker}
\title{Applying postblocking in a field}
\description{
  The blocker function takes a dataframe that contains columns indicating x and y coordinates and use them to create blocking with different numbers. For example, a square field with no blocks after applying the function will have several new columns for different blockings that are sequential, i.e. 2 blocks in one direction, then 4 blocks by blocking in the other direction and so on.
}
\usage{
  blocker(x,rows="ROWS",ranges="RANGES",by, char=FALSE)
}

\arguments{
    \item{x}{a dataframe with 2 obligatory columns; rows and ranges which can have different names and can be matched with the next arguments.}
  \item{rows}{the name of the numeric column that indicates one direction in the field.}
  \item{ranges}{the name of the numeric column that indicates the other direction in the field.}
  \item{by}{optional argument to indicate the name of the column of the dataframe x that indicates the environments so the field is filled by environment.}
  \item{char}{a TRUE/FALSE value indicating if the blocking column should be returned in a character class instead of numeric.}
  
}
\value{
  \describe{
    
    \item{$fin}{ a new dataframe with new columns for different blocking proposals.}
  }
}
\references{
  Fikret Isik. 2009. Analysis of Diallel Mating Designs. North Carolina State University, Raleigh, USA.
  
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
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
CPpheno <- CPdata$pheno#[,-c(1:4)]
CPgeno <- CPdata$geno
#### look at the data
head(CPpheno)
#### fill the design
gg <- fill.design(x=CPpheno, rows="Row",ranges="Col")
head(gg)
#### apply the postblocking 
gg2 <- blocker(x=gg, rows="Col",ranges="Row", char = TRUE)
head(gg2)
#### now you can use them in your mixed models

}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}