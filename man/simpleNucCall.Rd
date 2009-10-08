\name{simpleNucCall}
\alias{simpleNucCall}
%- Also NEED an '\alias' for EACH other topic documented here.
		\title{Predict nucleosome positions from high-throughput sequencing data}
\description{
	This function provides a simplified interface to \code{\link{callBindingSites}} with defaults suitable for
	the detection of nucleosomes.
}
\usage{
	simpleNucCall(data, bind=128, support=17, background=2000, chrLen, ...)
}
%- maybe also 'usage' for other objects documented here.
		\arguments{
			\item{data}{Either an object of class \code{\link[ShortRead:AlignedRead]{AlignedRead}} or a list. See below for details.}
			\item{bind}{Length of binding region to use (see Details).}
			\item{support}{Length of support region to use (see Details).}
			\item{background}{Length of background window. If this is missing it will be set to 10*(\code{bind}+2*\code{support}).}
			\item{chrLen}{Numeric vector indicating the length of all chromosomes. Only needed when 
				\code{data} is an \code{\link[ShortRead:AlignedRead]{AlignedRead}} object. \code{\link[ShortRead:readBfaToc]{readBfaToc}} 
				may be used to supply this information.}
			\item{...}{Further arguments to \code{\link{callBindingSites}}}
		}

\value{
	A \code{list} with components
	\item{binding}{A \code{data.frame} with columns \sQuote{chromosome}, \sQuote{position}, \sQuote{score} and \sQuote{pvalue}
		indicating the centre of predicted binding sites together with their score and associated p-value.}
	\item{score}{A \code{list} with all calculated scores. One numeric vector per chromosome.}
	\item{pval}{A \code{list} with all corrected p-values. One numeric vector per chromosome.}
}
\references{ ~put references to the literature/web site here ~ }
\author{Peter Humburg}

\seealso{\code{\link{callBindingSites}} for additional parameters.}
\examples{
## generate some simple artificial read data
set.seed(1)

## determine binding site locations
b <- sample(1:1e6, 5000)

## sample read locations
fwd <- unlist(lapply(b, function(x) sample((x-83):(x-73), 20, replace=TRUE)))
rev <- unlist(lapply(b, function(x) sample((x+73):(x+83), 20, replace=TRUE)))

## add some background noise
fwd <- c(fwd, sample(1:(1e6-25), 50000))
rev <- c(rev, sample(25:1e6, 50000))

## create data.frame with read positions as input to strandPileup
reads <- data.frame(chromosome="chr1", position=c(fwd, rev), 
	length=25, strand=factor(rep(c("+", "-"), times=c(150000, 150000))))

## create object of class ReadCounts
readPile <- strandPileup(reads, chrLen=1e6, extend=1, plot=FALSE)

## predict binding site locations
bindScore <- simpleNucCall(readPile, bind=147, support=20, plot=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
		% R documentation directory.
\keyword{models}
\keyword{htest}




