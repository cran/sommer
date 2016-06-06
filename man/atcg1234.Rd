\name{atcg1234}
\alias{atcg1234}
\title{
Letter to number converter
}
\description{
This function was designed to help users to transform their data in letter format to numeric format. Details in the format are not complex, just a dataframe with markers in columns and individuals in rows. Only markers, NO extra columns of plant names etc (names of plant can be stored as rownames).

}
\usage{
atcg1234(data, ploidy=2, format="ATCG", maf=0, multi=TRUE, silent=FALSE)
}
\arguments{
  \item{data}{a dataframe with markers in columns and individuals in rows. Preferable the rownames are the ID of the plants so you don't lose track of what is what.}
  \item{ploidy}{a numeric value indicating the ploidy level of the specie. The default is 2 which means diploid.}
  \item{format}{one of the two possible values allowed by the program "ATCG", which means your calls are in base-pair-letter code, i.e. "AT" in a diploid call, "AATT" tetraploid etc (just example). Therefore possible codes can be "A", "T", "C", "G", "-" (deletion), "+" (insertion). Alternatively "AB" format can be used as well. Commonly this depends from the genotyping technologies used, such as GBS or microarrays. In addition, we have enabled also the use of single-letter code used by Cornell, i.e. A=AA, C=CC, T=TT, G=GG, R=AG, Y=CT, S=CG, W=AT, K=GT, M=AC. The "ATCG" format also works for the bi-allelic marker codes from join map such as "lm", "ll","nn", "np","hh","hk","kk"}
  \item{maf}{minor allele frequency used to filter the SNP markers, the default is zero which means all markers are returned in numeric format.}
  \item{multi}{a TRUE/FALSE statement indicating if the function should get rid of the markers with more than 2 alleles. If FALSE, which indicates that if markers with multiple alleles are found, the alternate and reference alleles will be the first 2 alleles found. This could be risky since some alleles will be masked, i.e. AA AG AT would take only A and G as reference and alternate alleles, converting to numeric format 2 1 1, giving the same effect to AG and AT which could be a wrong assumption. The default is TRUE, removes markers with more than two alleles. }
  \item{silent}{a TRUE/FALSE value indicating if a progress bar should be drawn for each step of the conversion. The default is silent=FALSE, which means that we want progress bar to be drawn.}
  
}
\value{
\describe{

\item{$data}{ a new dataframe of markers in numeric format with markers in columns and individuals in rows.}

}
}

\author{
Giovanny Covarrubias-Pazaran
}

\examples{
data(PolyData)
genotypes <- PolyData$PGeno
genotypes[1:5,1:5] # look the original format

####=================================================####
#### convert markers to numeric format polyploid potatoes
####=================================================####
#numo <- atcg1234(data=genotypes, ploidy=4)
#numo[1:5,1:5]; dim(numo)

####=================================================####
#### convert markers to numeric format diploid rice lines
####=================================================####
#data(RICE)
#X <- RICE$RiceGeno; X[1:5,1:5]; dim(X)
#numo2 <- atcg1234(data=X, ploidy=2)
#numo2[1:5,1:5]; dim(numo2)
}
\seealso{The core function of the package \code{\link{mmer}}}