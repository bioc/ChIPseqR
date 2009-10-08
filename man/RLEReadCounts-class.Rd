\name{RLEReadCounts-class}
\Rdversion{1.1}
\docType{class}
\alias{RLEReadCounts-class}
\alias{chrLength,RLEReadCounts-method}
\alias{nreads,RLEReadCounts-method}

\title{Run-length Encoded Read Counts}
\description{This class provides a memory efficient representation of strand specific read counts.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{ReadCounts(counts, names, compress = TRUE)} or by calls
to \code{\link{strandPileup}}. 
}
\section{Slots}{
  \describe{
    \item{\code{counts}:}{Object of class \code{"list"} with one component per chromosome, 
    containing a read counts encoded in an object of class \code{\linkS4class{RleList}}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ReadCounts}"}, directly.
}
\section{Methods}{
  \describe{
    \item{chrLength}{\code{signature(x = "RLEReadCounts")}: Returns length of all chromosomes represented in \code{x}.}
%    \item{compress}{\code{signature(x = "RLEReadCounts")}: }
    \item{decompress}{\code{signature(x = "RLEReadCounts")}: Expands read counts and returns object of class \code{\linkS4class{ReadCounts}}.}
    \item{nreads}{\code{signature(x = "RLEReadCounts")}: Returns the number of reads on each chromosome, split by strand (if \code{byStrand} is \code{TRUE})}
    \item{plot}{\code{signature(x = "RLEReadCounts", y = "missing")}: Generates plots of read counts.}
	 }
}
\author{Peter Humburg}

\seealso{
	\code{\linkS4class{ReadCounts}}, \code{\linkS4class{RleList}}  
}
\examples{
showClass("RLEReadCounts")

## generate some very simple artificial read data
set.seed(1)
fwd <- sample(c(50:70, 250:270), 30, replace=TRUE)
rev <- sample(c(197:217, 347:417), 30, replace=TRUE)
## create data.frame with read positions as input to strandPileup
reads <- data.frame(chromosome="chr1", position=c(fwd, rev), 
	length=25, strand=factor(rep(c("+", "-"), times=c(30, 30))))

## create object of class ReadCounts
readPile <- strandPileup(reads, chrLen=500, extend=1, plot=FALSE, compress=TRUE)

names(readPile)
length(readPile)
sapply(readPile, sum)
}
\keyword{classes}
