\name{ReadCounts-class}
\Rdversion{1.1}
\docType{class}
\alias{ReadCounts-class}
\alias{[<-,ReadCounts,ANY,missing-method}
\alias{[,ReadCounts,ANY,missing,missing-method}
\alias{[[<-,ReadCounts,ANY,missing-method}
\alias{[[,ReadCounts,ANY,missing-method}
\alias{$<-,ReadCounts-method}
\alias{$<-,ReadCounts,ANY-method}
\alias{$,ReadCounts-method}
\alias{lapply,ReadCounts-method}
\alias{length<-,ReadCounts-method}
\alias{length,ReadCounts-method}
\alias{names<-,ReadCounts,ANY-method}
\alias{names,ReadCounts-method}
\alias{sapply,ReadCounts-method}
\alias{nreads,ReadCounts-method}
\alias{chrLength,ReadCounts-method}
\alias{ReadCounts}
\alias{nreads}
\alias{chrLength}

\title{Class "ReadCounts"}
\description{Represents counts of (possibly extended) reads for each strand of the genome.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{ReadCounts(counts, names, compress=FALSE)} or by calls
to \code{\link{strandPileup}}.
}
\section{Slots}{
	 \describe{
    \item{\code{counts}:}{Object of class \code{"list"} with one component per chromosome, containing a matrix of read counts (one column per strand).}
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "ReadCounts", i = "ANY", j = "missing")}: Replace read counts for chromosomes indicated by \code{i}. }
    \item{[}{\code{signature(x = "ReadCounts", i = "ANY", j = "missing", drop = "missing")}: Returns list of read counts for chromosomes indicated by \code{i}.}
    \item{[[<-}{\code{signature(x = "ReadCounts", i = "ANY", j = "missing")}: Replace read counts for chromosome \code{i}.}
    \item{[[}{\code{signature(x = "ReadCounts", i = "ANY", j = "missing")}: Returns read counts for chromosome \code{i}.}
    \item{$<-}{\code{signature(x = "ReadCounts")}: Replace read counts for chromosome \code{i} (by name).}
    \item{$}{\code{signature(x = "ReadCounts")}: Returns read counts for chromosome \code{i} (by name).}
    \item{callBindingSites}{\code{signature(data = "ReadCounts")}: Predict bindingsites from read counts.}
    \item{chrLength}{\code{signature(x = "ReadCounts", subset = "ANY")}: Returns length of all chromosomes represented in \code{x}.}
    \item{lapply}{\code{signature(X = "ReadCounts")}: Apply function to read counts for each chromosome.}
    \item{length<-}{\code{signature(x = "ReadCounts")}: Change the number of chromosomes represented by \code{x} to \code{value}.}
    \item{length}{\code{signature(x = "ReadCounts")}: Number of chromosomes represented by \code{x}.}
    \item{names<-}{\code{signature(x = "ReadCounts", value = "ANY")}: Change names of chromosomes.}
    \item{names}{\code{signature(x = "ReadCounts")}: Chromosome names.}
    \item{nreads}{\code{signature(x = "ReadCounts", byStrand = "Logical", subset = "ANY")}: Returns the number of reads on each chromosome, split by strand (if \code{byStrand} is \code{TRUE}).}
    \item{sapply}{\code{signature(X = "ReadCounts")}: Apply function to read counts for each chromosome.}
	 }
}
\usage{
ReadCounts(counts=list(), names=NULL, compress=TRUE)
}
\arguments{
\item{counts}{A list of read counts. Each component is a two column matrix of strand specific read counts for a chromosome.}
\item{names}{Character vector of chromosome names. If this is \code{NULL} the names of \code{counts} are used instead.}
\item{compress}{Logical indicating whether read counts should be compressed.}
}
\author{Peter Humburg}

\seealso{
	\code{\linkS4class{BindScore}}, \code{\link{strandPileup}} , \code{\link[ChIPseqR:compress-ReadCounts]{compress}}, \code{\link{decompress}}
}
\examples{
showClass("ReadCounts")

## generate some very simple artificial read data
set.seed(1)
fwd <- sample(c(50:70, 250:270), 30, replace=TRUE)
rev <- sample(c(197:217, 347:417), 30, replace=TRUE)
## create data.frame with read positions as input to strandPileup
reads <- data.frame(chromosome="chr1", position=c(fwd, rev), 
	length=25, strand=factor(rep(c("+", "-"), times=c(30, 30))))

## create object of class ReadCounts
readPile <- strandPileup(reads, chrLen=500, extend=1, plot=FALSE, compress=FALSE)

names(readPile)
length(readPile)
sapply(readPile, sum)

}
\keyword{classes}
