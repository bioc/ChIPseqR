\name{strandPileup}
\alias{strandPileup-methods}
\alias{strandPileup,AlignedRead-method}
\alias{strandPileup,data.frame-method}
\alias{strandPileup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Strand specific read counts}
\description{
  Given a set of aligned reads this function computes the number of reads
  starting at each position in the genome.
}
\usage{
\S4method{strandPileup}{AlignedRead}(aligned, chrLen, extend, compress = TRUE, plot = TRUE, ask = FALSE, ...)
\S4method{strandPileup}{data.frame}(aligned, chrLen, extend, compress = TRUE, plot = TRUE, ask = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aligned}{An object containing information about aligned reads (see Details).}
  \item{chrLen}{A numeric vector giving the length of each chromosome.}
  \item{extend}{A numeric value indicating how far reads should be extended.}
  \item{compress}{Logical indicating whether read counts should be compressed.}
  \item{plot}{If this is \code{TRUE} (the default) read coverage is plotted for all chromosomes.}
  \item{ask}{Logical. Setting this to \code{TRUE} causes the system to wait for user input before displaying a new plot.
  See \code{\link[grDevices:devAskNewPage]{devAskNewPage}}.}
  \item{\dots}{Further arguments to \code{\link[ShortRead:pileup]{pileup}}.}
}
\details{
  The method for \code{data.frame} requires the column names to follow a strict naming scheme. Required columns are
 \describe{
  \item{\sQuote{chromosome}}{A factor with chromosome names.}
  \item{\sQuote{strand}}{A factor with levels \dQuote{-} and \dQuote{+} indicating which strand the read mapped to.}
  \item{\sQuote{start} \bold{or} \sQuote{position}}{Start position of read on chromosome.}
  \item{\sQuote{end} \bold{or} \sQuote{length}}{End position of read or length of read respectively.}
  }
}

\value{
  An object of class \code{\linkS4class{ReadCounts}}.
}
\author{Peter Humburg}

\seealso{\code{\link[ShortRead:pileup]{pileup}}, \code{\link[ShortRead:AlignedRead]{AlignedRead}}, 
\code{\link{callBindingSites}}}
\examples{
	## generate some very simple artificial read data
	set.seed(1)
	fwd <- sample(c(50:70, 250:270), 30, replace=TRUE)
	rev <- sample(c(197:217, 347:417), 30, replace=TRUE)
	## create data.frame with read positions as input to strandPileup
	reads <- data.frame(chromosome="chr1", position=c(fwd, rev), 
		length=25, strand=factor(rep(c("+", "-"), times=c(30, 30))))
	
	readPile <- strandPileup(reads, chrLen=500, extend=1, plot=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
