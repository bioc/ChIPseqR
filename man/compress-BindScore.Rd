\name{compress-BindScore}
\alias{compress-BindScore}
\alias{compress,BindScore-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compress BindScore Objects
}
\description{
Generates a compressed representation of binding site scores.
}
\usage{
\S4method{compress}{BindScore}(x, digits=16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An object of class \code{\linkS4class{BindScore}}.
}
  \item{digits}{
	Integer indicating the number of decimal places to retain.
}
}
\details{
Binding site scores are compressed by first rounding them to \code{digits} decimal
places and then applying run-length encoding. 
}
\value{
An object of class \code{\linkS4class{RLEBindScore}}.
}

\author{
Peter Humburg
}
\note{
	Compression reduces the precision of binding site scores and may affect results, especially for small values
	of \code{digits}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{Rle}}, \code{\linkS4class{RleList}}, \code{\link{compress-BindScore}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
