\name{h2.fun}
\alias{h2.fun}
\title{
Obtain heritabilities with three different methods
}
\description{
Obtain heritabilities based on three different methods; Cullis et al. (2006), Oakey at al. (2006), and line-mean h2 (Falconer, 1995).

Cullis at al. (2006): h2= 1 - (PEM.mu/Vg)

Oakey at al. (2006): h2= 1 - (tr(G^-1 * Czz)/m)

Falconer (1995): h2 = Vg / (Vg+(Ve/r*e))

where "PEV.mu" is the average prediction error variance for the genetic term, "Vg" and "Ve" are the genetic and residual variance respectively estimated by REML, "G^-1" is the inverse of A*Vg where "A" is the additive relationship matrix, "Czz" is the prediction error variance for the genetic term, "m" is the number of test lines, "r" is the replicates per environment and "e" the number of environments.
}
\usage{

h2.fun(object, data, gTerm=NULL, eTerm=NULL)
      
}
\arguments{
  \item{object}{a model fitted with the mmer or mmer2 functions.}
  \item{data}{the dataset used to fit the model provided in the object argument.}
  \item{gTerm}{a character vector specifying the genetic terms fitted in the model.}
  \item{eTerm}{a character vector specifying the environment term fitted in the model.}
}
\details{
Please see the description or go to the canonical papers where methods are explained with more detail.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Oakey, Helena, et al. "Joint modeling of additive and non-additive genetic line effects in single field trials." Theoretical and Applied Genetics 113.5 (2006): 809-819.

Cullis, Brian R., Alison B. Smith, and Neil E. Coombes. "On the design of early generation variety trials with correlated data." Journal of Agricultural, Biological, and Environmental Statistics 11.4 (2006): 381-393.

Falconer, Douglas S., Trudy FC Mackay, and Richard Frankham. "Introduction to quantitative genetics (4th edn)." Trends in Genetics 12.7 (1996): 280.
}
\examples{

data(h2example)
head(h2example)
####=========================================####
#### fit the mixed model (very heavy model)
####=========================================####
# ans1 <- mmer2(y~Env,
#               random=~at(Env):Name + at(Env):Block,
#               rcov=~at(Env):units,
#               data=h2example)
# summary(ans1)

####=========================================####
#### calculate heritability
####=========================================####
# length(levels(h2example$Env)) # number of environment
# 2 # number of blocks
# # using the pin function
# pin(ans1,h2example~(V1+(V3/15))/(V1+(V3/15)+(V5/(2*15))))
# # using the 
# h2.fun(ans1, h2example, gTerm = "Name", eTerm = "Env")


}
\seealso{
\code{\link{sommer}}, \code{\link{mmer2}}
}