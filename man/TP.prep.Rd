\name{TP.prep}
\alias{TP.prep}
\title{Selecting the best training population for genomic selection}
\description{
This function was designed to help users build the best training population for real genomic selection applications. The \code{\link{RICE}} dataset from Zhao et al. (2011) has been incorporated to the package to highlight this new feature of sommer implemented in v.1.5. The TP.prep function is just a function that uses genetic markers or any kind of numeric matrix to perform principal component analysis, identify the names of the validation population the user provide, and find the most suitable training population for performing genomic selection based on similarities or similarities-dissimilarities.

My recomendation would be to randomly choose the TP and predict the performance of the VP. Do this a 100 times and get the average GEBV for the VP and take decisions using those.

}
\usage{
TP.prep(markers=NULL, vp.names=NULL, tp.size=seq(50,200,25), method="sim-dissim", 
        npop=300, nelite=1, mutprob=.5, niterations=100, lambda=NULL)

}
\arguments{
  \item{markers}{a marker matrix in numeric format but any numeric matrix can be used in this argument. No extra columns should be in this matrix. ROWNAMES SHOOULD BE THE NAMES OF THE INDIVIDUALSIN THE POPULATION.}
  \item{vp.names}{a character vector with the names of the validaton population. The marker matrix should have rownames with the individuals so the program can find them and identify the most suitable training population.}
  \item{tp.size}{a numeric vector indicating the population sizes desired for the training population.}
  \item{method}{a method among "similar", "sim-dissim", "random", "PEV", or "CDmean" indicating if the program should look for similar plants (based on the PCA of the marker matrix) or a mixture of similar and constrasting individuals, select plants at random, minimize the PEV, or calculate the CD mean. The default is "sim-dissim", meaning will make a mixture.}
  \item{npop}{genetic algorithm parameter, number of solutions at each iteration.}
  \item{nelite}{genetic algorithm parameter, number of solutions selected as elite parents which will generate the next set of solutions.}
  \item{mutprob}{genetic algorithm parameter, probability of mutation for each generated solution.}
  \item{niterations}{genetic algorithm parameter, number of iterations.}
  \item{lambda}{a value for the CDmean algorithm where lambda=Var(e)/Var(g).}

}

\value{
If all parameters are correctly indicated the program will return:
\describe{

\item{$tp.list}{a list with the names of the TP populations. Each element of the list corresponds to each population size specified in the 'tp.size' argument.}

}
}
\references{
Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds, Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in Oryza sativa. Nat Comm 2:467 DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011.

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
data(RICE)
W <- RICE$RiceGenoN; W[1:5,1:5]; dim(W)
Y <- RICE$RicePheno; Y[1:5,1:5]; dim(Y)
(vp.str <- RICE$vp.str[1:25]) # validation population with structure
bases <- seq(50,200,25) # TP sizes required

####=========================================####
#### a) Random sample based
####=========================================####
rest <- (1:dim(W)[1])[-which(rownames(W) \%in\% vp.str)]
#tp.plant.strR <- TP.prep(markers=W, vp.names = vp.str, 
#                         tp.size = bases, method="random")

####=========================================####
#### b) Genetic simmilarity-dissimilarity based
####=========================================####
#tp.plant.strG1 <- TP.prep(markers=W, vp.names = vp.str, 
#                         tp.size = bases, method="sim-dissim")
                         
#tp.plant.strG2 <- TP.prep(markers=W, vp.names = vp.str, 
#                         tp.size = bases, method="sim")
                       
####=========================================####
#### c) QTL based
#### we assume real situation that VP 
#### wouldn't be phenotyped for GWAS
####=========================================####
#K1 <- A.mat(W)
#Z1 <- diag(dim(K1)[1])
#image(K1)
#ETA1 <- list(list(Z=Z1,K=K1))
#YY <- Y; YY[vp.str,] <- NA
#ans.GWAS <- mmer(Y=YY$Protein.content, Z=ETA1, W=W)
#(h2 <- sqrt(ans.GWAS$var.comp[1,]/sum(ans.GWAS$var.comp)))
#X <- bag(ans.GWAS, nmar = 25); head(X); dim(X)

#tp.plant.strG3 <- TP.prep(markers=X, vp.names = vp.str, 
#                      tp.size = bases, method="sim-dissim")

####=========================================####
#### COMPARE METHODS with 50 plants
####=========================================####
#metho <- list(a=tp.plant.strR[[1]], b=tp.plant.strG1[[1]],
#      c=tp.plant.strG2[[1]], d=tp.plant.strG3[[1]])

#resos <- numeric()
#for(i in 1:4){ # for each method
  
#  tpgi <- metho[[i]] # training pop genetic-based
#  Wpgi <- W[c(vp.str,tpgi),] # provisional markers based on genes
#  Ypgi <- Y[c(vp.str,tpgi),] # provisional phenos based on genes
#  Xpgi <- X[c(vp.str,tpgi),] # provisional bag-mat based on genes
  
#  Kpgi <- A.mat(Wpgi)
#  Zpgi <- diag(dim(Kpgi)[1])
#  Ypgi[vp.str,] <- NA
#  yt <- Ypgi$Protein.content
#  ETAg <- list(list(Z=Zpgi, K=Kpgi))
#  mixg <- mmer(Y=yt, X=Xpgi, Z=ETAg, method="EMMA") #X=Xpgi,
#  resos[i] <- cor(Y[vp.str,"Protein.content"], 
#                   fitted(mixg)[(1:length(vp.str)),], 
#                   use="complete")
#}
#resos
### the random method needs more iterations to have a real idea what
### would be the predictive ability of using plants at random,
### this is just to give an idea to users.
}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}