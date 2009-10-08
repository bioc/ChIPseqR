\name{plot,ReadCounts,missing-method}
\alias{plot,ReadCounts,missing-method}
\alias{plot,RLEReadCounts,missing-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic Plots for Read Counts
}
\description{
Produces plots to assess the distribution of reads, either for an entire chromosome or within
a (small) window.
}
\usage{
\S4method{plot}{ReadCounts,missing}(x, chr, center, score, width=2000, type=c("hilbert", "window"), ...)
\S4method{plot}{RLEReadCounts,missing}(x, chr, center, score, width=2000, type=c("hilbert", "window"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	Object of class \code{\linkS4class{ReadCounts}} or \code{\linkS4class{RLEReadCounts}}. 
}
  \item{chr}{
	Index or name of chromosome for which read counts should be plotted.
}
  \item{center}{
	For type \sQuote{window}, the center coordinate of the window to plot.  
}
  \item{score}{
  	For type \sQuote{window}, an object of type \code{\linkS4class{BindScore}} 
  	(or \code{\linkS4class{BindScore}}) that should be used to include information
  	about the score and predicted binding sites into the plot.
}
  \item{width}{
	For type \sQuote{window}, the width of the window.  
}
  \item{type}{
	Character string indicating the type of plot that should be produced (see details).  
}
  \item{\dots}{
	Further arguments (see details).
}
}
\details{
Type \sQuote{window} produces a plot that shows read counts as vertical bars. If \code{score}
is not missing, it is used to plot the score and predicted binding sites (if any) as well. 
Any additional arguments are passed on to \code{\link{.plotWindow}}.

Type \sQuote{hilbert} produces a Hilbert curve plot of read counts for an entire chromosome.
Additional arguments are passed on to  \code{\link{.plotReads}}.
}
\value{
Called for its side effect.
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[HilbertVis:hilbertImage]{hilbertImage}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
