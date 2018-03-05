\name{overlay}
\alias{overlay}
\title{Overlay Matrix}
\description{
  `overlay` adds r times the design matrix for model term t to the existing design matrix. Specifically, if the model up to this point has p effects and t has a effects, the a columns of the design matrix for t are multiplied by the scalar r (default value 1.0). This can be used to force a correlation of 1 between two terms as in a diallel analysis.
}
\usage{
  overlay(..., rlist=NULL, prefix=NULL)
}

\arguments{
  \item{...}{as many vectors to factor vectors to do the overlay.}
  \item{rlist}{a list of scalar values indicating the times that each incidence matrix overlayed should be multiplied by. By default r=1.}
  \item{prefix}{a character name to be added before the column names of the final overlay matrix. This may be useful if you have entries with names starting with numbers which programs such as asreml doesn't like, or for posterior extraction of parameters, that way 'grep'ing is easier.}
}
\value{
  \describe{
    
    \item{$S3}{ an incidence matrix with as many columns levels in the vectors provided to build the incidence matrix.}
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
modh <- mmer2(sugar~1, random=~overlay(female,male) + geno, 
              data=HDdata)
summary(modh)
#### model using overlay and covariance structures

# A <- diag(7); A[1,2] <- 0.5; A[2,1] <- 0.5 # fake covariance structure
# colnames(A) <- as.character(1:7); rownames(A) <- colnames(A);A
# 
# modh2 <- mmer2(sugar~1, random=~ overlay(female,male) + geno,
#                G=list(`overlay(female, male)`=A),data=HDdata)
# summary(modh2)
# 
# ####=========================================####
# #### Example 2 (matrix-based)
# #### using overlay with mmer function
# ####=========================================####
# 
# #### GCA matrix for half diallel using male and female columns
# #### use the 'overlay' function to create the half diallel matrix
# Z1 <- overlay(HDdata$male, HDdata$female)
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