\name{atcg1234}
\alias{atcg1234}
\title{
Letter to number converter
}
\description{
This function was designed to help users to transform their data in letter format to numeric format. Details in the format are not complex, just a dataframe with markers in columns and individuals in rows. Only markers, NO extra columns of plant names etc (names of plant can be stored as rownames).

}
\usage{
atcg1234(data, ploidy=2, format="ATCG", maf=0)
}
\arguments{
  \item{data}{a dataframe with markers in columns and individuals in rows. Preferable the rownames are the ID of the plants so you don't lose track of what is what.}
  \item{ploidy}{a numeric value indicating the ploidy level of the specie. The default is 2 which means diploid.}
  \item{format}{one of the two possible values allowed by the program "ATCG", which means your calls are in base-pair-letter code, i.e. "AT" in a diploid call, "AATT" tetraploid etc (just example). Alternatively "AB" format can be used as well. Commonly this depends from the genotyping technologies used, such as GBS or microarrays. In addition, we have enabled also the use of single-letter code used by Cornell, i.e. A=AA, C=CC, T=TT, G=GG, R=AG, Y=CT, S=CG, W=AT, K=GT, M=AC.}
  \item{maf}{minor allele frequency used to filter the SNP markers, the default is zero which means all markers are returned in numeric format.}
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
####### convert markers to numeric format
#numo <- atcg1234(data=genotypes, ploidy=4)
#numo[1:5,1:5]; dim(numo)
}