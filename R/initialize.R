# constructors for S4 classes defined in classes.R
# 
# Author: Peter Humburg
###############################################################################

setMethod("initialize", c("ReadCounts"),
		function(.Object, counts=list(), names=NULL, ...){
			stopifnot(is.list(counts))
			
			.Object <- callNextMethod(.Object, ...)
			
			## ensure we have chromosome names
			if(!is.null(names)) names(counts) <- names
			if(is.null(names(counts))) 
				stop("Chromosome names are missing.")
			
			## ensure that counts is integer matrix
			if(!all(sapply(counts, is.matrix))){
				counts <- lapply(counts, 
						function(x){
							if(is.matrix(x)) return(x)
							else{
								warning("Converting entry to two column matrix.")
								return(matrix(x, ncol=2))
							} 
						})
			}
			if(!all(sapply(counts, is.integer))){
				counts <- lapply(counts, function(x) matrix(as.integer(x), ncol=2))
			}
			
			## set column names to indicate strand
			counts <- lapply(counts, "colnames<-", c("+", "-"))
	
			.Object@counts <- counts
			
			.Object
		}
)

setMethod("initialize", c("RLEReadCounts"), 
		function(.Object, counts=list(), names=NULL, ...){
			stopifnot(is.list(counts))
			## ensure we have chromosome names
			if(!is.null(names)) names(counts) <- names
			if(is.null(names(counts))) 
				stop("Chromosome names are missing.")
			
			## ensure that counts is integer matrix
			idx <- !(sapply(counts, is.matrix) | sapply(counts, is, "RleList"))
			if(any(idx)){
				counts[idx] <- lapply(counts[idx], 
						function(x){
							warning("Converting entry to two column matrix.")
							return(matrix(x, ncol=2)) 
						})
			}
			idx <- !(sapply(counts, is.integer) | sapply(counts, is, "RleList"))
			if(any(idx))
				counts[idx] <- lapply(counts[idx], function(x) matrix(as.integer(x), ncol=2))
			
			## compress read counts
			idx <- !sapply(counts, is, "RleList")
			counts[idx] <- lapply(counts[idx], function(x) IRanges::RleList(IRanges::Rle(x[, 1]), IRanges::Rle(x[, 2])))
			counts[idx] <- lapply(counts[idx], "names<-", c("+", "-"))
			
			.Object@counts <- counts
			
			.Object
		}
)

## assumes that the function that was called to compute these results had arguments
## 'bind' and 'support' referring to the length of binding and support regions
setMethod("initialize", "BindScore",
		function(.Object, functionCall, score=list(), pvalue=list(), peaks=list(), cutoff=c(-Inf, 1), 
				nullDist=c(0, 1), names=NULL, start=1L, digits=16, ...){
			.Object <- callNextMethod(.Object, ...)
			
			## ensure we have chromosome names
			if(any(sapply(list(names(score), names(pvalue), names(peaks)), is.null)) ||
					!(all.equal(names(score), names(pvalue)) && all.equal(names(score), names(peaks)))){
				## ensure all names match given chromosome names
				names(score) <- names
				names(pvalue) <- names
				names(peaks) <- names
			}
			if(is.null(names(score))) 
				stop("Chromosome names are missing.")
			
			if(length(cutoff) < 2) stop("Argument 'cutoff' has to be of length 2.")
			if(length(cutoff) > 2){
				warning("'cutoff' has length ", length(cutoff), ". Truncating to length 2.")
				length(cutoff) <- 2
			}
			names(cutoff) <- c("score", "pvalue")
			
			if(length(nullDist) < 2) stop("Argument 'nullDist' has to be of length 2.")
			if(length(nullDist) > 2){
				warning("'nullDist' has length ", length(nullDist), ". Truncating to length 2.")
				length(nullDist) <- 2
			}
			names(nullDist) <- c("mean", "sd")
			
			.Object@functionCall <- functionCall
			.Object@score <- score
			.Object@pvalue <- pvalue
			.Object@peaks <- peaks
			.Object@cutoff <- cutoff
			.Object@nullDist <- nullDist
			.Object@start <- as.integer(start)
			
			.Object
		}
)

setMethod("initialize", "RLEBindScore",
		function(.Object, functionCall, score=list(), pvalue=list(), peaks=list(), cutoff=c(-Inf, 1), 
				nullDist=c(0, 1), names=NULL, start=1L, digits=16, ...){
			
			## ensure we have chromosome names
			if(any(sapply(list(names(score), names(pvalue), names(peaks)), is.null)) ||
					!(all.equal(names(score), names(pvalue)) && all.equal(names(score), names(peaks)))){
				## ensure all names match given chromosome names
				names(score) <- names
				names(pvalue) <- names
				names(peaks) <- names
			}
			if(is.null(names(score))) 
				stop("Chromosome names are missing.")
			
			if(length(cutoff) < 2) stop("Argument 'cutoff' has to be of length 2.")
			if(length(cutoff) > 2){
				warning("'cutoff' has length ", length(cutoff), ". Truncating to length 2.")
				length(cutoff) <- 2
			}
			names(cutoff) <- c("score", "pvalue")
			
			if(length(nullDist) < 2) stop("Argument 'nullDist' has to be of length 2.")
			if(length(nullDist) > 2){
				warning("'nullDist' has length ", length(nullDist), ". Truncating to length 2.")
				length(nullDist) <- 2
			}
			names(nullDist) <- c("mean", "sd")
			
			.Object@functionCall <- functionCall
			
			## compress scores and p-values
			idx <- !sapply(score, is, "Rle")
			score[idx] <- lapply(score[idx], function(x){ 
						x[is.na(x)] <- -Inf 
						x <- Rle(round(x, digits=digits))
						runValue(x)[is.infinite(runValue(x))] <- NA
						x
					}
			) 
			.Object@score <- score
			
			idx <- !sapply(pvalue, is, "Rle")
			pvalue[idx] <- lapply(pvalue[idx], function(x){ 
						x[is.na(x)] <- Inf 
						x <- Rle(round(x, digits=16))
						runValue(x)[is.infinite(runValue(x))] <- NA
						x
					}
			)
			.Object@pvalue <- pvalue  
						
			.Object@peaks <- peaks
			.Object@cutoff <- cutoff
			.Object@nullDist <- nullDist
			.Object@start <- as.integer(start)
			
			.Object
		}
)


## constructors that allow to choose between compressed and non-compressed versions
## These could also contain some logic to determine whether compression is likely to be beneficial
ReadCounts <- function(counts=list(), names=NULL, compress=TRUE){
	if(compress) new("RLEReadCounts", counts=counts, names=names)
	else new("ReadCounts", counts=counts, names=names)
}

BindScore <- function(call, score=list(), pvalue=list(), peaks=list(), cutoff=c(-Inf, 1), 
		nullDist=c(0, 1), names=NULL, start=1L, compress=TRUE, digits=16){
	if(compress) new("RLEBindScore", functionCall=call, score=score, pvalue=pvalue, peaks=peaks, cutoff=cutoff,
				nullDist=nullDist, names=names, start=start, digits=digits)
	else new("BindScore", functionCall=call, score=score, pvalue=pvalue, peaks=peaks, cutoff=cutoff,
				nullDist=nullDist, start=start, names=names)
}
