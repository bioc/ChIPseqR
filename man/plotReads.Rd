\name{plotReads}
\alias{.plotReads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot compact representation of read counts on a chromosome}
\description{
  Creates an image of all read counts for a chromosome.
}
\usage{
.plotReads(x, scale = c("total", "ratio"), log = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A two column matrix with read counts (as produced by \code{\link{strandPileup}}).}
  \item{scale}{Character string indicating whether to plot the total number of reads on both strands or the ratio between them.}
  \item{log}{Logical indicating whether to use log read counts.}
  \item{\dots}{Further plotting arguments.}
}
\details{
 The read counts (or read count ratios) are plotted as a Hilbert curve using 
 \code{\link[HilbertVis:hilbertImage]{hilbertImage}}.
}
\value{
  Called for its side effect.
}
\references{\cite{Anders, Simon, 'Visualization of genomic data with the Hilbert curve',  Bioinformatics , vol. 25, no. 10, 1231-1235 (2009).}}
\author{Peter Humburg}

\seealso{\code{\link{strandPileup}}, \code{\link[HilbertVis:hilbertImage]{hilbertImage}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{internal}
