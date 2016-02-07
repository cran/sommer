\name{hdm}
\alias{hdm}
\title{Half Diallel Matrix}
\description{
This function uses a dataframe with 2 columns named "female" and "male" in numeric format and creates an incidence matrix for a single explanatory variable corresponding to the GCA effect. The resulting incidence matrix can be used in the mmer function as a 'Z' argument in the 2-level list argument for random effects.
}
\usage{
hdm(data)
}

\arguments{
  \item{data}{a dataframe with 2 columns named 'female' and 'male' with numeric or factor values indicating the male or female used to produce such hybrid.}
}
\value{
\describe{

\item{$Z}{ an incidence matrix with as many columns as parents in the dataframe indicating with ones the parents used for a particular hybrid (in rows).}
}
}
\references{
Fikret Isik. 2009. Analysis of Diallel Mating Designs. North Carolina State University, Raleigh, USA.
}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
data(HDdata)
head(HDdata)
#### GCA matrix for half diallel using male and female columns
#### use the 'hdm' function to create the half diallel matrix
#Z1 <- hdm(HDdata[,c(3:4)])
#### SCA matrix
#Z2 <- model.matrix(~as.factor(geno)-1, data=HDdata)
#### response variable
#y <- HDdata$sugar
#ETA <- list(list(Z=Z1), list(Z=Z2)) # Zu component
#modHD <- mmer(y=y, Z=ETA, draw=FALSE, silent=TRUE)
#summary(modHD)
}