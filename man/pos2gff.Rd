\name{pos2gff}
\alias{pos2gff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert genome coordinates into GFF format}
\description{
  Provides facility to export the location of genomic features to a GFF formatted file.
}
\usage{
pos2gff(pos, method, feature, len, strand, score, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{Named list with one component per chromosome giving the start position of features on that chromosome.}
  \item{method}{Entry for method field in GFF file. Recycled as necessary}
  \item{feature}{Entry for feature field in GFF file. Recycled as necessary}
  \item{len}{Length of fetures. This is used to calculate matching end positions for each start position given in \code{pos}}
  \item{strand}{Entry for feature field in GFF file. Recycled as necessary}
  \item{score}{Entry for feature field in GFF file. Recycled as necessary}
  \item{name}{Entry for feature field in GFF file. Recycled as necessary}
}
\value{
  A \code{data.frame} with columns \code{'chromosome'}, \code{'method'}, \code{'feature'}, \code{'start'}, 
  \code{'end'}, \code{'score'}, \code{'strand'}. Writing this data frame to a text 
  file produces a GFF formatted file.  
}
\references{The GFF specification: http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml}
\author{Peter Humburg}
\examples{
pos <- list(chr1=c(10, 50, 60), chr2=c(22, 200, 500))
pos2gff(pos, "test", "foo", 25, c("+", "+", "-", "+", "-", "-"), 0, "test")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
