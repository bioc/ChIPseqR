\name{getCutoff}
\alias{getCutoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine significance threshold for binding site scores}
\description{
  Given a vector of observed binding site scores and a desired false discovery rate, this
  function returns the lowest score that should be considered significant to achieve the 
  given false discovery rate.
}
\usage{
getCutoff(score, alpha = 0.05, tailCut = 0.95, adapt = FALSE, lambda, plot = TRUE, returnPval = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{score}{Numeric vector with binding site scores.}
  \item{alpha}{Desired false discovery rate.}
  \item{tailCut}{Truncation point used to exclude outliers when fitting the null distribution.}
  \item{adapt}{Logical indicating whether an adaptive false discovery rate should be used.}
  \item{lambda}{If \code{adapt} is \code{TRUE} this is used in estimating the proportion of 
  scores that is unrelated to binding sites.}
  \item{plot}{If this is \code{TRUE} (the default) a plot of the observed score distribution
  and estimated null distribution is generated.}
  \item{returnPval}{Indicates whether or not the corrected p-values for all scores should be returned.}
}

\value{
  A list with components
  \item{cutoff}{A numeric vector with the score and nominal false discovery rate corresponding
  to the determined cutoff.}
  \item{h0}{A numeric vector with the mean and standard deviation of the estimated null distribution.}
  \item{pvalue}{If \code{returnPval} is \code{TRUE}, the p-values corresponding to the scores in \code{score}. 
  Note that all missing values are removed.}
  \item{pi0}{If \code{adapt} is \code{TRUE}, the estimated proportion of scores not related to binding sites.}
}

\note{
This function is used by \code{\link{callBindingSites}} to determine the significance threshold for
peak-calling.
}

\references{For the adaptive false discovery rate procedure used if \code{adapt=TRUE} see 
JD Storey, JE Taylor and D Siegmund. Strong control, conservative point
estimation and simultaneous conservative consistency of false discovery rates: a unified approach.
Journal of the Royal Statistical Society B, 66(1):187-205, 2004.
}

\author{Peter Humburg}

\seealso{\code{\link{callBindingSites}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
