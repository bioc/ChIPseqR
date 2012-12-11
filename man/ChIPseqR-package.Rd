\name{ChIPseqR-package}
\alias{ChIPseqR-package}
\alias{ChIPseqR}
\docType{package}
\title{
Identifying Protein Binding Sites in High-Throughput Sequencing Data
}
\description{
ChIPseqR provides a set of functions for the analysis of ChIP-seq data. Protein binding sites
are located by identifying a characteristic pattern of peaks in read counts on both DNA strands.  
}
\details{
\tabular{ll}{
Package: \tab ChIPseqR\cr
Type: \tab Package\cr
Version: \tab 1.13.1\cr
Date: \tab 2012-12-11\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

The easiest way to obtain binding site predictions for nucleosomes is to use \code{\link{simpleNucCall}}. This provides
a simple interface to \code{\link{callBindingSites}}. This function operates on \code{\link[ShortRead:AlignedRead]{AlignedRead}} 
objects and provides useful defaults for nucleosome analysis. Parameters can be adjusted to detect the presence of other 
DNA binding proteins, e.g. transcription factors. If more fine control is desired the following steps will produce binding site 
predictions:
\describe{
\item{\code{\link{strandPileup}}:}{Turn mapped reads into read counts along the genome.}
\item{\code{\link{startScore}}:}{Score potential binding sites.}
\item{\code{\link{getCutoff}}:}{Determine cutoff required to achieve desired false discovery rate.}
\item{\code{\link{pickPeak}}:}{Find all peaks in the binding site score that exceed the significance 
	threshold determined by \code{getCutoff}. These are the predicted binding sites.}  
}
}
\author{
Peter Humburg

Maintainer: Peter Humburg <Peter.Humburg@well.ox.ac.uk>
}
\references{
Humburg, P., Helliwell, C., Bulger, D. & Stone, G. ChIPseqR: Analysis of ChIP-seq Experiments. BMC Bioinformatics 12, 39+ (2011).
}

\keyword{ package }
\seealso{
\code{\link[ShortRead:ShortRead-package]{ShortRead}}
}
\examples{
	## See 'simpleNucCall' for examples of how to obtain nucleosome predictions. 
}
