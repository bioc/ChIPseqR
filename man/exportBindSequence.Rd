\name{exportBindSequence}
\Rdversion{1.1}
\alias{exportBindSequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Export sequence of predicted binding sites
}
\description{
Extracts sequence of predicted binding sites from reference genome and exports them in FASTA format.
}
\usage{
exportBindSequence(prediction, reference, bind, overlap = FALSE, file = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prediction}{
Object of class \code{\linkS4class{BindScore}}.
}
  \item{reference}{
Reference genome sequence (as \code{\linkS4class{XStringSet}} object). 
}
  \item{bind}{
Length of binding site to assume for sequence extraction. This may be missing in which case the value is 
derived from \code{'prediction'}.
}
  \item{overlap}{
Logical indicating whether overlapping predictions should be allowed.
}
  \item{file}{
Name of output file.
}
}

\value{
An \code{\linkS4class{XStringViews}} object containing the sequences. If a file name is provided this is returned
invisibly.
}
\references{
Package \code{Biostrings}
}
\author{
Peter Humburg
}


\seealso{
\code{\linkS4class{XStringViews}}, \code{\linkS4class{XStringSet}}, \code{\linkS4class{BindScore}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
