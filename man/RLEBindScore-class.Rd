\name{RLEBindScore-class}
\Rdversion{1.1}
\docType{class}
\alias{RLEBindScore-class}


\title{Run-length Encoded Binding Site Scores}
\description{This class provides a memory efficient representation of binding site scores.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{BindScore(functionCall, score, pvalue, peaks, cutoff, nullDist, names, start, digits, compress=TRUE)}
or through calls to \code{\link{callBindingSites}}. 
}
\section{Slots}{
  \describe{
    \item{\code{functionCall}:}{Object of class \code{"call"} storing the function call used to initiate the analysis.}
    \item{\code{score}:}{Object of class \code{"list"}. The binding site score. One run-length encoded numeric vector per chromosome.}
    \item{\code{pvalue}:}{Object of class \code{"list"}. The (adjusted and run-length encoded) p-values corresponding to the scores in slot \code{score}.}
    \item{\code{peaks}:}{Object of class \code{"list"} giving the location of significant peaks in the binding site score. These correspond to the location of predicted binding sites.}
    \item{\code{cutoff}:}{Object of class \code{"numeric"} with entries \sQuote{pvalue} and \sQuote{score} giving the significance threshold used for peak calling in terms of p-value and score.}
    \item{\code{nullDist}:}{Object of class \code{"numeric"} providing the parameters of the null distribution used to determine p-values.}
    \item{\code{start}:}{Object of class \code{"integer"} indicating the index corresponding to the first entry in \code{score} (assumed to be the same for all chromosomes).}
  }
}
\section{Extends}{
Class \code{"\linkS4class{BindScore}"}, directly.
}
\section{Methods}{
  \describe{
    %\item{compress}{\code{signature(x = "RLEBindScore")}: ... }
    \item{decompress}{\code{signature(x = "RLEBindScore")}: conversion to \code{\linkS4class{BindScore}} object.}
	 }
}
\author{Peter Humburg}

\seealso{
\code{\linkS4class{BindScore}}, \code{\linkS4class{Rle}}
}
\examples{
showClass("RLEBindScore")

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
bindScore <- simpleNucCall(readPile, bind=147, support=20, plot=FALSE, compress=TRUE)

## number of binding sites found
length(bindScore)

## the first few predictions, by score
head(bindScore)

## score and p-value cut-off used
cutoff(bindScore)

}
\keyword{classes}
