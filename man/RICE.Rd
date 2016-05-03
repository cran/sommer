\name{RICE}
\alias{RICE}
\title{wheat lines dataset}
\description{

Information from a collection of 413 rice lines.  The RICE data set is from 
Rice Diversity Org. Program. The lines are genotyped with 36,901 SNP markers and phenotyped for more than 30 traits. This data set was included in the package to illustrate the selection of training populations for obtainin maximum prediction accuracies. We have found that the best way to select a training population is either 1) taking a random sample from the entire population, 2) do PCA and selecting very similar and very different genotypes (quite similar to random sampling),  or 3) performing a GWAS with the data available to identify regions of significance and do PCA but using only markers that are in significant regions, (selecting from plants close to the VP and also far, to make sure we sample from different haplotypes for the QTLs). All versions are shown in the examples with rice data using the \code{\link{TP.prep}} function.

}

\usage{
  data(RICE)
}

\format{
 RicePheno contains the phenotypes
 RiceGeno contains genotypes letter code
 RiceGenoN contains the genotypes in numerical code using atcg1234 converter function
}

\source{
  Rice Diversity Organization
  http://www.ricediversity.org/data/index.cfm.
}

\references{

Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds, Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in Oryza sativa. Nat Comm 2:467 DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011.

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
#rest <- (1:dim(W)[1])[-which(rownames(W) \%in\% vp.str)]
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
#ans.GWAS <- mmer(y=YY$Protein.content, Z=ETA1, W=W)
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
#  mixg <- mmer(y=yt, X=Xpgi, Z=ETAg, method="EMMA") #X=Xpgi,
#  resos[i] <- cor(Y[vp.str,"Protein.content"], 
#                   fitted(mixg)[(1:length(vp.str)),], 
#                   use="complete")
#}
#resos
### the random method needs more iterations to have a real idea what
### would be the predictive ability of using plants at random,
### this is just to give an idea to users.
}
\keyword{datasets}
