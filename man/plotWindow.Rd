\name{plotWindow}
\alias{.plotWindow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot read counts within a genomic region}
\description{
  Read count within a selected region of the genome are plotted, optionally together 
  with binding site score and location of predicted binding sites.
}
\usage{
.plotWindow(data, chr, center, score, width=2000, bind, start, end, bind.col=3, score.type='l',	
	xlab=NULL, ylab="Read count", cutoff=TRUE, offset = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Object of class \code{\linkS4class{ReadCounts}} or a list of read counts.}
  \item{chr}{A character string or numeric index identifying the the chromosome on which the region is located.}
  \item{center}{Numeric value giving the center of the region on chromosome \code{'chr'} that should be plotted.}
  \item{score}{An object of class \code{\linkS4class{BindScore}}, may be missing.}
  \item{width}{Width of the window to plot. The plotted region will be [center - width/2, center + width/2].}
  \item{bind}{Length of binding site, ignored if \code{'score'} is missing.}
  \item{start}{Start of plotting window (may be used together with \code{end} instead of \code{center}).}
  \item{end}{End of plotting window (may be used together with \code{start} instead of \code{center}).}
  \item{bind.col}{Color used to indicate location of binding sites, ignored if \code{'score'} is missing.}
  \item{score.type}{Plotting type to use for score, ignored if \code{'score'} is missing.}
  \item{xlab}{X-axis label. This defaults to a description of the genomic location constructed from \code{'chr'},
  \code{'center'} and \code{'width'}}
  \item{ylab}{Y-axis label.}
  \item{cutoff}{Logical indicating whether the significance threshold used to predict binding site locations 
  should be indicated by a horizontal line in the plot.}
  \item{offset}{Position on chromosome \code{chr} corresponding to the first read count in \code{data}.}
  \item{\dots}{Further arguments to \code{plot}.}
}
\details{
  If \code{'score'} is present the binding site score is plotted on top of read counts. Scores are rescaled
  to lie between 0 and the maximum number of reads in the window. A corresponding scale is added to the right 
  hand axis.
}
\value{
  Called for its side effect.
}
\author{Peter Humburg}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{internal}