\name{acessors}
\Rdversion{1.1}
\alias{peaks}
\alias{peaks-methods}                              
\alias{peaks,ANY-method}
\alias{pvalue}
\alias{pvalue-methods}                              
\alias{pvalue,ANY-method}
\alias{score}
\alias{cutoff}
\alias{cutoff-methods}                              
\alias{cutoff,ANY-method}
\alias{cutoff<-}
\alias{cutoff<--methods}                              
\alias{cutoff<-,ANY-method}
\alias{support}
\alias{support-methods}                              
\alias{support,ANY-method}
\alias{binding}
\alias{binding-methods}                              
\alias{binding,ANY-method}
\alias{nullDist}
\alias{nullDist-methods}                              
\alias{nullDist,ANY-method}
\alias{nullDist<-}
\alias{nullDist<--methods}                              
\alias{nullDist<-,ANY-method}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Access slots of S4 classes
}
\description{
Accessor functions for S4 classes in package \code{"ChIPseqR"}. 
}
\usage{
	\S4method{binding}{ANY}(x, ...)
	\S4method{cutoff}{ANY}(x, ...)
	\S4method{cutoff}{ANY}(x, ...) <- value
	\S4method{nullDist}{ANY}(x, ...)
	\S4method{nullDist}{ANY}(x, ...) <- value
	\S4method{peaks}{ANY}(x, ...)
	\S4method{pvalue}{ANY}(x, ...)
	\S4method{support}{ANY}(x, ...)
}
\section{Methods}{
\describe{                  
\item{x = "ANY"}{Default method for accessor function.}
}}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An S4 object.
}
  \item{\dots}{
Further arguments, ignored by default method.
}
\item{value}{New value for slot.}
}
\details{
These methods allow safe read (and in some cases write) access to slots of S4 classes and should be used
for this purpose rather than modifying slots manually.
}
\value{
The current value of the interrogated slot.
}

\author{
Peter Humburg
}
\note{
This page documents the generics and their default behaviour. See the help page of each class for class 
specific implementations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{ReadCounts}}, \code{\linkS4class{BindScore}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

