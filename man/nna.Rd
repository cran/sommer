\name{nna}
\alias{nna}

\title{
Nearest neighbour adjustment
}
\description{
This function takes a dataframe that contains the variables 'row', 'col', and 'y' to create a new variable 'nnx' which is a new covariable to be used in the linear model based on the neighboring rows and columns. By default it uses the formula used by Lado et al. (2013) where the new variable is defined as:

nnx_i = yi - mean(y1,y2,...,y6)

where spatially such plants will look like:

[...] [....] [y2] [...] [...]

[y5] [y1] [yi] [y3] [y6]

[...] [....] [y4] [...] [...]

Corners and edge plots are treated by only using the neighbor plots available, i.e.:

[y2] [...] [...]

[yi] [y3] [y4]

[y1] [...] [...]

The number of plants to use in the adjustment can be modified with the additional arguments 'nrows' and 'ncols'. Once each new 'y' is computed a new variable named 'nnx' is added to the original dataframe which can be used as an extra covariate in the model.

}
\usage{
nna(pheno,trait="y",rown="row",coln="col",nrows=1,ncols=2)
}

\arguments{
  \item{pheno}{a dataframe that contains the variables 'row', 'col', and 'y'.}
  \item{trait}{an additional name for the response variable in case the user wants to avoid to change the names of the dataframe.}
  \item{rown}{an additional name for the row variable in case the user wants to avoid to change the names of the dataframe.}
  \item{coln}{an additional name for the column variable in case the user wants to avoid to change the names of the dataframe.}
  \item{nrows}{the number of row neighbors to be used in the adjustment for the plants.}
  \item{ncols}{the number of column neighbors to be used in the adjustment for the plants.}
}
\details{

The function assumes unreplicated data, one measure per plot and a single experiment. If the user has a dataframe with the same data replicated in different environments or blocks this should be separated so the function do not adjust for neighbouring plots in other years, repetitions or environments.

}
\value{
If everything is defined correctly the function returns:
\describe{
\item{$pheno}{the original dataframe with a new variable named 'nnx' which contains the nearest neighbor adjustment.}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Lado et al. (2013) Increased genomic prediction accuracy in wheat breeding through spatial adjustment of field trial data. G3, 3:2105:2114

}

\examples{

data(yates.oats)

head(yates.oats)
newyates <- nna(yates.oats, trait="Y")
head(newyates)

plot(newyates$Y, newyates$nnx)
cor(newyates$Y, newyates$nnx)

#### now fit the models and compare #####

m3 <- mmer2(fixed=Y ~ V*N, random = ~ B + B:MP,
            data = yates.oats)
yates.oats$res <- residuals(m3)

m4 <- mmer2(fixed=Y ~ V*N + nnx, random = ~ B + B:MP,
            data = newyates)
newyates$res <- residuals(m4)

# library(lattice)
# wireframe(res~row*col,yates.oats)
# wireframe(res~row*col,newyates)

}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmer2}}}
