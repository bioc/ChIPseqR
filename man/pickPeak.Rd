\name{pickPeak}
\alias{pickPeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify peaks above a given threshold}
\description{
  Given a vector of scores and a threshold, this function finds all peaks 
  that exceed the threshold.
}
\usage{
pickPeak(score, threshold, offset = 0, sub = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{score}{Numeric vector.}
  \item{threshold}{All values in \code{score} below this value are ignored.}
  \item{offset}{Offset to add to the determined peak locations.}
  \item{sub}{Logical. If this is \code{FALSE} (the default) for each region that exceeds the threshold
  only the global maximum is returned. Otherwise local maxima are returned as well.}
}

\value{
  If \code{sub = FALSE} a numeric vector giving the location of all peaks.
  Otherwise a list with components
  \item{peaks}{The same peak locations that are returned for \code{sub = FALSE}.}
  \item{subPeaks}{A list with one component for each entry in \sQuote{peaks} giving the location
  of local maxima.} 
}

\note{
This function is used by \code{\link{callBindingSites}} for peak-calling. 
}
\author{Peter Humburg}

\seealso{\code{\link{callBindingSites}}, \code{\link{startScore}}, \code{\link{getCutoff}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
