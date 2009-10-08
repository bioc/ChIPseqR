\name{plot-BindScore}
\alias{plot,BindScore,missing-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic Plots for Binding Site Scores 
}
\description{
Generates plots to assess the fit of the estimated null distribution.
}
\usage{
\S4method{plot}{BindScore,missing}(x, npoints = 10000, type=c("density", "qqplot"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An object of class \code{\linkS4class{BindScore}}.
}
  \item{npoints}{
	Maximum number of points to plot in a QQ-plot.
}
   \item{type}{
	Character string indicating the plot type.  
}
  \item{\dots}{
 Further arguments (currently ignored).
}
}
\details{
Type \sQuote{density} produces a histogram of binding site scores with overlaid null distribution.
Type \sQuote{qqplot} produces a normal QQ-plot comparing the observed binding site scores to the 
null distribution.
}
\value{
Called for its side effect.
}

\author{
Peter Humburg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
