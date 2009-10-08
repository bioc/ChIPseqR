\name{windowCounts}
\alias{windowCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize read counts in a sliding window}
\description{
  Read counts are summarized in a sliding window of variable size with variable overlap between windows.
}
\usage{
windowCounts(reads, window = 1000, shift = 500, method = sum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{Numeric vector of read counts.}
  \item{window}{Width of window.}
  \item{shift}{Distance between consecutive window start positions.}
  \item{method}{Function used to produce a summary for each window. It should accept a single numeric vector as argument.}
}

\value{
  If \code{method} returns a single value a vector of all window summaries is returned, otherwise
  the return value is a list with one component for each window.
}
\author{Peter Humburg}

\examples{
## generate some very simple artificial read data
set.seed(1)
fwd <- sample(c(50:70, 250:270), 30, replace=TRUE)
rev <- sample(c(197:217, 347:417), 30, replace=TRUE)
## create data.frame with read positions as input to strandPileup
reads <- data.frame(chromosome="chr1", position=c(fwd, rev), 
	length=25, strand=factor(rep(c("+", "-"), times=c(30, 30))))

## create object of class ReadCounts
readPile <- strandPileup(reads, chrLen=501, extend=1, plot=FALSE, compress=FALSE)

## get number of reads in sliding window
wdwCount <- windowCounts(apply(readPile[[1]], 1, sum), window=10, shift=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
