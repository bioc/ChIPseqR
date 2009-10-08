\name{decompress}
\alias{decompress}
\alias{decompress,RLEBindScore-method}
\alias{decompress,RLEReadCounts-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Read Count and Binding Site Score Representations 
}
\description{
These methods extract read count and binding site sores from compressed representations. 
}
\usage{
\S4method{decompress}{RLEReadCounts}(x)
\S4method{decompress}{RLEBindScore}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An object of class \code{\linkS4class{RLEBindScore}} or \code{\linkS4class{RLEReadCounts}}.
}
}

\value{
An object of class \code{\linkS4class{BindScore}} or \code{\linkS4class{ReadCounts}} respectively.
}
\author{
Peter Humburg
}

\seealso{
\code{\link[=compress-methods]{compress}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
