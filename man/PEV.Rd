\name{PEV}
\alias{PEV}
\title{Selecting the best training population for genomic selection}
\description{
This function is a wrapper from the STPGA package to obtain the best subset of individuals to predict a group of individuals minimizing the predicted error variance (PEV). Is used internally in the TP.prep function.

}
\usage{
PEV(PCAs,candidates,Test,ntoselect, npop, nelite, mutprob, niterations, lambda)

}
\arguments{
  \item{PCAs}{nxk PCA matrix from the predictor variables.}
  \item{candidates}{vector of names for the population without the test set.}
  \item{Test}{name of the individuals in the test or VP set.}
  \item{ntoselect}{number of individuals to select in the training population.}
  \item{npop}{number of solutions at each iteration.}
  \item{nelite}{number of elite solutions for TP pops.}
  \item{mutprob}{probability of mutation for each solution}
  \item{niterations}{number of iterations}
  \item{lambda}{scalar shrinkage in PEV.}

}

\value{
If all parameters are correctly indicated the program will return:
\describe{

\item{$tp.list}{a list with the names of the TP populations. Each element of the list corresponds to each population size specified in the 'tp.size' argument.}

}
}
\references{
Akdemir, Deniz. "Training population selection for breeding value" prediction. 2014.

}
\author{
Giovanny Covarrubias-Pazaran
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}