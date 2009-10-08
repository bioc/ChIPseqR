\name{callBindingSites-methods}
\docType{methods}
\alias{callBindingSites-methods}
\alias{callBindingSites}
\alias{callBindingSites,ANY-method}
\alias{callBindingSites,character-method}
\alias{callBindingSites,matrix-method}
\alias{callBindingSites,ReadCounts-method}
\title{Predict protein binding sites from high-throughput sequencing data}
\description{
 Methods for function \code{callBindingSites} in Package `ChIPseqR'. These methods are used to 
 identify protein binding sites from ChIP-seq data. 
}
\section{Methods}{
\describe{

\item{data = "ANY"}{Default method to handle all forms of input not explicitly handled by their own method. In particular this will be used for objects of class \code{\link[ShortRead:AlignedRead]{AlignedRead}} and \code{data.frame} but it will handle class for which a \code{\link{strandPileup}} method is available.}

\item{data = "character"}{Allows to use a file name referring to a file of mapped sequence reads as input.}

\item{data = "matrix"}{Uses a matrix of read counts (for a single chromosome) as input.}

\item{data = "ReadCounts"}{This methods implements the peak calling algorithm. Other methods will typically reformat their input and pass it on to this method.}
}}

\usage{
\S4method{callBindingSites}{ANY}(data, chrLen, plot=TRUE, verbose=TRUE, ..., plotTo)
\S4method{callBindingSites}{character}(data, type, minQual=70, ...)
\S4method{callBindingSites}{matrix}(data, chrName="chr", ...)
\S4method{callBindingSites}{ReadCounts}(data, bind, support, background, bgCutoff=0.9, supCutoff=0.9, 
fdr = 0.05, extend=1, tailCut=0.95, piLambda=0.5, adapt=FALSE, corSummary=median, compress = TRUE,
digits = 16, plot=TRUE, verbose=TRUE, ask=FALSE, plotTo, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Either an object containing information about mapped reads or a list. See below for details.}
  \item{bind}{Length of binding region to use (see Details).}
  \item{support}{Length of support region to use (see Details).}
  \item{background}{Length of background window. If this is missing it will be set to 10*(\code{bind}+2*\code{support}).}
  \item{chrLen}{Numeric vector indicating the length of all chromosomes. Only needed when 
  \code{data} is an \code{\link[ShortRead:AlignedRead]{AlignedRead}} object. \code{\link[ShortRead:readBfaToc]{readBfaToc}} 
  may be used to supply this information.}
  \item{bgCutoff}{Numeric value between 0.5 and 1. This determines how much estimates of the background read density are 
  allowed to vary for adjacent windows. Set to 1 to disable cutoff.}
  \item{supCutoff}{Numeric value between 0.5 and 1. This determines how much estimates of the support region read density are 
  allowed to vary for forward and reverse strand. Set to 1 to disable cutoff.}
  \item{fdr}{Target false discovery rate.}
  \item{extend}{Numeric value indicating how far mapped reads should be extended when calculating read counts.}
  \item{type}{Format of alignment file (see \code{\link[ShortRead:readAligned]{readAligned}} forr details).}
  \item{minQual}{Minimum alignment quality to use. All reads with lower alignment quality are discarded.}
  \item{tailCut}{Truncation point used to exclude outliers when estimating null distribution.}
  \item{chrName}{Name to use for the single chromosome.}
  \item{piLambda}{If \code{adapt=TRUE} this parameter is used to estimate the proportion of scores not related to binding sites.}
  \item{adapt}{Logical indicating whether an adaptive false discovery rate should be used. If this is \code{FALSE} (the default)
  the usual Benjamini-Hochberg procedure is used to control the FDR.}
  \item{corSummary}{Function used to summarise cross-correlation across chromosomes. See the Details section on 
  binding and support region. }
  \item{compress}{Logical indicating whether the return value should be compressed.}
  \item{digits}{Number of decimal places to retain for binding site score for compression.}
  \item{plot}{Logical. If \code{plot=TRUE} (the default) some diagnostic plots are produced during the analysis.}
  \item{verbose}{Logical. If \code{verbose=TRUE} (the default) status messages are printed to indicate progress.}
  \item{ask}{Logical. Setting this to \code{TRUE} causes the system to wait for user input before displaying a new plot.
  See \code{\link[grDevices:devAskNewPage]{devAskNewPage}}.}
  \item{plotTo}{Character string giving the name of a file that should be used to store plots generated during the analysis.
  If this is not missing a pdf file with the given name will be created.}
  \item{\dots}{Additional arguments. Most methods pass them on to the \code{ReadCounts} method.}
}
\details{
  The length of binding and support regions can either be given as a single value or as a range of 
  possible values (by providing the minimum and maximum). In the latter case the cross-correlation
  between read counts on forward and reverse strand will be used to determine a value within that range. 
  Note that this may lead  sub-optimal choices of binding and support region length.
}
\value{
  An object of class \code{\linkS4class{BindScore}} if \code{compress = FALSE}, otherwise an object
  of class \code{\linkS4class{RLEBindScore}}
}
\seealso{\code{\link{simpleNucCall}} for an interface with nucleosome specific defaults. This function uses 
\code{\link{strandPileup}}, \code{\link{startScore}}, \code{\link{getCutoff}} and \code{\link{pickPeak}}. See the 
help pages of these functions for additional detail on the individual steps involved. See \code{\link{getBindLen}} 
for details on the estimation of binding and support region length.}

\examples{
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
bindScore <- callBindingSites(readPile, bind=147, support=20, background=2000, plot=FALSE)
}

\keyword{methods}
\keyword{models}
\keyword{htest}
