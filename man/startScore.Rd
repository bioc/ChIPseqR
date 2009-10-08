\name{startScore}
\alias{startScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Score potential protein binding sites}
\description{
  For each position in the genome this function computes a score indicating the
  likelihood that a protein binding site starts at that position.
}
\usage{
startScore(data, b, support, background, bgCutoff, supCutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A two column matrix with read counts. The two columns correspond to reads on the 
  forward and reverse strand respectively.}
  \item{b}{Length of binding region.}
  \item{support}{Length of support region.}
  \item{background}{Length of background window.}
  \item{bgCutoff}{Cutoff for the change in read rates between adjacent windows (see Details).}
  \item{supCutoff}{Cutoff for the change in read rates between support regions on forward and 
  reverse strand (see Details).}
}
\details{
  Robust estimates of read rates in background windows and support regions are obtained by limiting the
  difference between related estimates. Consider a forward support region of length 10 containing 20 reads.
  The maximum likelihood estimate for the rate parameter of the (assumed) underlying Poisson distribution is
  \eqn{\hat{\lambda} = \frac{20}{10} = 2}{lambda_hat = 20/10 = 2}. If there are 50 reads in the reverse
  support region a robust estimate of the rate parameter is obtained as 
  \code{max(50/10, qpois(supCutoff, lambda=lambda_hat))}  
}
\value{
  Numeric vector with binding site scores.
}
\note{
Instead of calling this function directly use \code{\link{callBindingSites}}.
}

\references{}
\author{Peter Humburg}

\seealso{\code{\link{callBindingSites}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
