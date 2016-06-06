\name{maxi.qtl}
\alias{maxi.qtl}
\title{Peak search by first derivatives}

\description{
This function find maximum peaks in a data frame from QTL mapping.
}
\usage{
maxi.qtl(sma, distan=10, q95=2.5)
}

\arguments{
 \item{sma}{A data frame with 3 basic columns; "chr", "pos", "lod". Make sure your data frame has this 3 column names. Even though not all the times the scale is measure as LOD this will not affect the peak search, is just to help the program recognize which column has the signal intensity measurement.}

\item{distan}{A numeric value indicating how far in distance (e.g. cM) we will consider a close QTL a different QTL since several points can form a QTL peak. The column "pos" is used to find such distances.}

\item{q95}{A numeric value indicating the threshold where we want to declare a peak significant in your scale (e.g. LOD, -log10(p.val), etc.).}

}
\details{
No major details. Is just a maximum peak searcher controlling for distance.
}
\value{
Retuns the biggest peaks for a data frame with chromosomes, position and intensities.
\describe{
\item{out}{ a data frame with the positions of the highest peaks.}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}
\seealso{The core function of the package: \code{\link{mmer}}}

