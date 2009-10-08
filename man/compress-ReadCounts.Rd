\name{compress-ReadCounts}
\alias{compress,ReadCounts-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compress ReadCount Objects
}
\description{
Generates a compressed representation of read counts. 
}
\usage{
\S4method{compress}{ReadCounts}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An object of class \code{\linkS4class{ReadCounts}}
}

}
\details{
Run-length encoding is used to obtain a compressed representation of read counts. 
}
\value{
An object of class \code{\linkS4class{RLEReadCounts}}
}

\author{
Peter Humburg
}

\seealso{
\code{\linkS4class{Rle}}, \code{\linkS4class{RleList}}, \code{\link{compress-BindScore}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

