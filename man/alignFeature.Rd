\name{alignFeature}
\alias{alignFeature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read counts relative to annotated features}
\description{
  Creates a set of (strand specific) read counts centred at the genomic features provided. 
}
\usage{
alignFeature(data, anno, offset = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{List with read counts as returned by \code{\link{strandPileup}}.}
  \item{anno}{Data frame with annotation data in GFF format.}
  \item{offset}{Half width of window around start point of annotated features.}
}
\value{
  List with one component for each feature in \code{anno}. 
}
\references{The GFF file format specification: http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml}
\author{Peter Humburg}

\examples{
set.seed(1)

## determine binding site locations
b <- sample(1:8.5e5, 500)

## sample read locations
fwd <- unlist(lapply(b, function(x) sample((x-83):(x-73), 20, replace=TRUE)))
rev <- unlist(lapply(b, function(x) sample((x+73):(x+83), 20, replace=TRUE)))

## add some background noise
fwd <- c(fwd, sample(1:(1e6-25), 5000))
rev <- c(rev, sample(25:1e6, 5000))

## create data.frame with read positions as input to strandPileup
reads <- data.frame(chromosome="chr1", position=c(fwd, rev), 
	length=25, strand=factor(rep(c("+", "-"), times=c(15000, 15000))))

## create object of class ReadCounts
readPile <- strandPileup(reads, chrLen=1e6, extend=1, plot=FALSE)

## convert binding site locations into GFF format
gff <- data.frame(chromosome="chr1", source="test", feature="binding", start=b-73, end=b+73, 
	score=".", strand=".", frame=".")

## align read counts relative to binding site location
aligned <- alignFeature(readPile, gff, offset=500)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
