\name{getBindLen}
\alias{getBindLen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate length of binding and support region}
\description{
  The cross-correlation between forward and reverse strand read counts is used to estimate the distance between
  peaks on both strands. This is then used to derive suitable values for the length of binding and support regions.
}
\usage{
getBindLen(data, bind, support, summary = median, verbose = FALSE, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An object of class \code{\linkS4class{ReadCounts}}.}
  \item{bind}{Either known length of binding region or minimum and maximum of binding region length to consider.}
  \item{support}{Either known length of support region or minimum and maximum of support region length to consider.}
  \item{summary}{Function to use to summarise cross-correlation across chromosomes.}
  \item{verbose}{Logical indicating whether progress messages should be printed.}
  \item{plot}{Logical indicating whether cross-correlation should be plotted.}
  \item{\dots}{Further arguments to \code{\link{getBindCor}}.}
}
\details{
 This assumes that the first peak in cross-correlation corresponds to the length of the binding site.
 Note that this is not accurate. The peak is closer to bind + 2*m where m is the median of the 
 read distribution in the support region ('read distribution in the support region' means the read 
 density as a function of distance to binding site start/end). Consequently this method will 
 overestimate the length of the binding site. If either bind or support are of length 1 
 this is assumed to be the known value and a more accurate estimate for the remaining parameter is used.
}
\value{
  A numeric vector giving the estimated lengths of binding and support regions.
}
\author{Peter Humburg}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

