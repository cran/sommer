\name{hdm}
\alias{hdm}
\title{Half Diallel Matrix}
\description{

THIS FUNCTION HAS BEEN OVERPASSED BY THE 'overlay' FUNCTION WHICH IS MORE FLEXIBLE AND FASTER.

This function uses a dataframe with 2 columns named "female" and "male" in numeric format and creates an incidence matrix for a single explanatory variable corresponding to the GCA effect. The resulting incidence matrix can be used in the mmer function as a 'Z' argument in the 2-level list argument for random effects.
}
\usage{
hdm(dat)
}

\arguments{
  \item{dat}{a dataframe with 2 columns named 'female' and 'male' with numeric or factor values indicating the male or female used to produce such hybrid.}
}
\value{
\describe{

\item{$Z}{ an incidence matrix with as many columns as parents in the dataframe indicating with ones the parents used for a particular hybrid (in rows).}
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
#### remove them and run the examples
####=========================================####
data(HDdata)
head(HDdata)

####=========================================####
#### Example 1 (formula-based)
#### using overlay with mmer2 function
####=========================================####
data(HDdata)
head(HDdata)
HDdata$female <- as.factor(HDdata$female)
HDdata$male <- as.factor(HDdata$male)
HDdata$geno <- as.factor(HDdata$geno)
#### model using overlay
modh <- mmer2(sugar~1, random=~female + and(male) + geno, 
              data=HDdata)
summary(modh)
#### model using overlay [and(.)] and covariance structures [g(.)]

# A <- diag(7); A[1,2] <- 0.5; A[2,1] <- 0.5 # fake covariance structure
# colnames(A) <- as.character(1:7); rownames(A) <- colnames(A);A
# 
# modh2 <- mmer2(sugar~1, random=~g(female) + and(g(male)) + geno, 
#                G=list(female=A, male=A),data=HDdata)
# summary(modh2)
# 
# ####=========================================####
# #### Example 2 (matrix-based)
# #### using overlay with mmer function
# ####=========================================####
# 
# #### GCA matrix for half diallel using male and female columns
# #### use the 'overlay' function to create the half diallel matrix
# Z1 <- overlay(HDdata[,c(3:4)])
# 
# #### Obtain the SCA matrix
# Z2 <- model.matrix(~as.factor(geno)-1, data=HDdata)
# 
# #### Define the response variable and run
# y <- HDdata$sugar
# ETA <- list(list(Z=Z1), list(Z=Z2)) # Zu component
# modHD <- mmer(Y=y, Z=ETA, draw=FALSE, silent=TRUE)
# summary(modHD)

}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}