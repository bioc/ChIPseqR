\name{BindScore-class}
\Rdversion{1.1}
\docType{class}
\alias{BindScore-class}
\alias{as.data.frame,BindScore-method}
\alias{[,BindScore,ANY,missing,missing-method}
\alias{[[,BindScore,ANY,missing-method}
\alias{[[,BindScore,ANY,numeric-method}
\alias{binding,BindScore-method}
\alias{cutoff<-,BindScore-method}
\alias{cutoff,BindScore-method}
\alias{head,BindScore-method}
\alias{lapply,BindScore-method}
\alias{length<-,BindScore-method}
\alias{length,BindScore-method}
\alias{max,BindScore-method}
\alias{min,BindScore-method}
\alias{names<-,BindScore,ANY-method}
\alias{names,BindScore-method}
\alias{nullDist,BindScore-method}
\alias{nullDist<-,BindScore-method}
\alias{peaks,BindScore-method}
\alias{range,BindScore-method}
\alias{score,BindScore-method}
\alias{support,BindScore-method}
\alias{tail,BindScore-method}
\alias{chrLength,BindScore-method}      
\alias{BindScore}               


\title{Class "BindScore"}
\description{This class provides the infrastructure to store results of ChIP-seq analysis.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BindScore", functionCall, score, pvalue, peaks, cutoff, nullDist, names, ...)}.
	 Objects of this class are typically created (and returned) by functions that perform
	 peak calling on ChIP-seq data. Usually there should be no need to create them directly. 
}
\section{Slots}{
	 \describe{
    \item{\code{functionCall}:}{Object of class \code{"call"} storing the function call used to initiate the analysis.}
    \item{\code{score}:}{Object of class \code{"list"}. The binding site score. One numeric vector per chromosome.}
    \item{\code{pvalue}:}{Object of class \code{"list"}. The (adjusted) p-values corresponding to the scores in slot \code{score}.}
    \item{\code{peaks}:}{Object of class \code{"list"} giving the location of significant peaks in the binding site score. These correspond to the location of predicted binding sites.}
    \item{\code{cutoff}:}{Object of class \code{"numeric"} with entries \sQuote{pvalue} and \sQuote{score} giving the significance threshold used for peak calling in terms of p-value and score.}
    \item{\code{nullDist}:}{Object of class \code{"numeric"} providing the parameters of the null distribution used to determine p-values.}
    \item{\code{start}:}{Object of class \code{"integer"} indicating the index corresponding to the first entry in \code{score} (assumed to be the same for all chromosomes).}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "BindScore")}: Convert results into a \code{data.frame} giving the location, score and p-value of significant peaks.}
    \item{[}{\code{signature(x = "BindScore", i = "ANY", j = "missing", drop = "missing")}: Restrict results to a subset of chromosomes. Chromosomes can either be identified by name or by numerical index.}
    \item{[[}{\code{signature(x = "BindScore", i = "ANY", j = "missing")}: Restrict results to a single chromosome. Note that \code{x[["chr1"]]} is identical to \code{x["chr1"]}.}
    \item{[[}{\code{signature(x = "BindScore", i = "ANY", j = "numeric")}: subset results to restrict them to a region on a single chromosome.}
    \item{binding}{\code{signature(x = "BindScore")}: Returns length of binding site used during analysis.}
    \item{chrLength}{\code{signature(x = "BindScore", subset = "ANY")}: Returns length of all chromosomes represented in \code{x}.}
    \item{cutoff<-}{\code{signature(x = "BindScore")}: Sets the significance cut-off. Argument \code{type=c("score", "pvalue")} determines which cut-off is to be set, the other is adjusted accordingly. This recalculates the significance of peaks in the binding site score and may be slow.}
    \item{cutoff}{\code{signature(x = "BindScore")}: Returns significance threshold used for analysis.}
    \item{head}{\code{signature(x = "BindScore")}: Returns the first \code{n} peaks. Argument \code{by = c("score", "position")} determines whether results are sorted by score or by genomic location.}
    \item{lapply}{\code{signature(X = "BindScore")}: Applies a function to results for each chromosome.}
    \item{length<-}{\code{signature(x = "BindScore")}: Reduces the number of chromosomes for which results are stored, i.e., \code{length(x) <- 3} only retains the first three chromosomes.}
    \item{length}{\code{signature(x = "BindScore")}: Returns the number of binding sites identified by the analysis.}
    \item{max}{\code{signature(x = "BindScore")}: Returns maximum score.}
    \item{min}{\code{signature(x = "BindScore")}: Returns minimum score.}
    \item{names<-}{\code{signature(x = "BindScore", value = "ANY")}: Sets the chromosome names.}
    \item{names}{\code{signature(x = "BindScore")}: Returns the chromosome names.}
    \item{nullDist<-}{\code{signature(x = "BindScore")}: Sets the parameters of the null distribution adjusting the significance cut-off in the process and predicts binding sites using the new null distribution.}
    \item{peaks}{\code{signature(x = "BindScore")}: Returns list of predicted binding sites.}
    \item{range}{\code{signature(x = "BindScore")}: Range of binding site scores.}
    \item{score}{\code{signature(x = "BindScore")}: Returns list of binding site scores.}
    \item{support}{\code{signature(x = "BindScore")}: Returns length of support region used during analysis.}
    \item{tail}{\code{signature(x = "BindScore")}: Returns the last \code{n} peaks. Argument \code{by = c("score", "position")} determines whether results are sorted by score or by genomic location.}
	 }
}
\usage{
	\S4method{binding}{BindScore}(x)
	\S4method{chrLength}{BindScore}(x, subset)
	\S4method{cutoff}{BindScore}(x, type=c("score", "pvalue"))
	\S4method{cutoff}{BindScore}(x, type=c("score", "pvalue")) <- value
	\S4method{head}{BindScore}(x, n=6, by=c("score", "position"), ...)
	\S4method{lapply}{BindScore}(X, FUN, ...)
	\S4method{length}{BindScore}(x)
	\S4method{length}{BindScore}(x) <- value
	\S4method{max}{BindScore}(x, ..., na.rm=TRUE)
	\S4method{min}{BindScore}(x, ..., na.rm=TRUE)
	\S4method{names}{BindScore}(x)
	\S4method{names}{BindScore,ANY}(x) <- value
	\S4method{nullDist}{BindScore}(x)
	\S4method{nullDist}{BindScore}(x) <- value
	\S4method{peaks}{BindScore}(x, ...)
	\S4method{range}{BindScore}(x, ..., na.rm=TRUE)
	\S4method{score}{BindScore}(x)
	\S4method{support}{BindScore}(x)
	\S4method{tail}{BindScore}(x, n=6, by=c("score", "position"), ...)
	BindScore(call, score=list(), pvalue=list(), peaks=list(), cutoff=c(-Inf, 1), nullDist=c(0, 1), names=NULL, start=1L, compress=TRUE, digits=16)
}
\arguments{
	\item{x}{Object of class \code{BindScore}.}
	\item{X}{Object of class \code{BindScore}.}
	\item{subset}{Index vector indicating a subset of \code{x}. If \code{subset} is  missing everything is selected.}
	\item{type}{A string indicating which type of cut-off should be returned or changed. Either \code{"score"} or \code{"pvalue"}}
	\item{n}{Number of entries to show.}
	\item{by}{A string indicating whether the output should be sorted by score or by position in the genome.}
	\item{na.rm}{Logical indication whether \code{NA}s should be ignored.}
	\item{FUN}{Function to apply to results for each chromosome.}
	\item{value}{Replacement value.}
	\item{call}{Function call used to generate the values of the other slots.}
	\item{score}{List of binding site scores. One component per chromosome.}
	\item{pvalue}{List of binding site p-values. One component per chromosome.}
	\item{peaks}{List of significant peaks in binding site score. One component per chromosome.}
	\item{cutoff}{Numeric vector of length two indicating the significance cut-off in terms of score and p-value.}
	\item{nullDist}{Parameters of the null distribution.}
	\item{names}{Character vector providing names for chromosomes.}
	\item{start}{Integer indicating position of first binding site score.}
	\item{compress}{Logical indicating whether scores and p-values should be compressed.}
	\item{digits}{The number of decimal places to retain for compression.}
	\item{\dots}{Further arguments.}
}
\references{ ~put references to the literature/web site here ~ }
\author{Peter Humburg}

\seealso{
	\code{\linkS4class{ReadCounts}} for the data structure used as input for the analysis and \code{\link[=callBindingSites,ReadCounts-method]{callBindingSites}}
}
\examples{
showClass("BindScore")

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
## the artificial dataset is very small so predictions may not be very reliable
bindScore <- simpleNucCall(readPile, bind=147, support=20, plot=FALSE, compress=FALSE)

## number of binding sites found
length(bindScore)

## the first few predictions, by score
head(bindScore)

## score and p-value cut-off used
cutoff(bindScore)

}
\keyword{classes}
