\name{decompress-methods}
\docType{methods}
\alias{decompress-methods}
\alias{decompress,ANY-method}
\alias{decompress,Rle-method}
\alias{decompress,RleList-method}

\title{Methods for Function decompress in Package `ChIPseqR'}
\description{
 Methods to extract compressed data structures.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{The default method simply returns \code{x}.}

\item{\code{signature(x = "Rle")}}{Restores the atomic vector encoded in \code{x}.}

\item{\code{signature(x = "RLEBindScore")}}{Returns an object of \code{\linkS4class{BindScore}}.}

\item{\code{signature(x = "RleList")}}{Extracts the components of \code{x} and restores them to atomic vectors.}

\item{\code{signature(x = "RLEReadCounts")}}{Returns an object of \code{\linkS4class{ReadCounts}}.}
}}
\keyword{methods}
