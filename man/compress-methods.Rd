\name{compress-methods}
\docType{methods}
\alias{compress}
\alias{compress-methods}
\alias{compress,RLEBindScore-method}
\alias{compress,RLEReadCounts-method}
\title{Methods for Function compress in Package `ChIPseqR'}
\description{
 Methods to obtain compressed versions of data structures. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "BindScore")}}{Converts \code{x} into an object of class \code{\linkS4class{RLEBindScore}}.}

\item{\code{signature(x = "ReadCounts")}}{Converts \code{x} into an object of class \code{\linkS4class{RLEReadCounts}}.}

\item{\code{signature(x = "RLEBindScore")}}{Returns (the already compressed) \code{x}.}

\item{\code{signature(x = "RLEReadCounts")}}{Returns (the already compressed) \code{x}.}
}}
\keyword{methods}
