\name{getBindCor}
\alias{getBindCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate cross-correlation between read counts}
\description{
  This function calculates the cross-correlation between read counts on forward and reverse strand.
}
\usage{
getBindCor(data, max.lag, summary, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An object of class \code{\linkS4class{ReadCounts}}}
  \item{max.lag}{Maximum lag to use in cross-correlation calculation.}
  \item{summary}{Function to use to summarise cross-correlation across chromosomes.}
  \item{plot}{Logical indicating whether to plot cross-correlation.}
  \item{\dots}{Further arguments, currently ignored.}
}
\details{
 Function \code{\link[timsac:fftcor]{fttcor}} in package \dQuote{timsac} is used to carry out the computation.
}
\value{
 The (summarised) cross-correlation. If \code{summary} is missing a list of cross-correlations for
 each chromosome is returned.
}
\author{Peter Humburg}

\seealso{\code{\link[timsac:fftcor]{fftcor}}, \code{\linkS4class{ReadCounts}}, \code{\link{getBindLen}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

