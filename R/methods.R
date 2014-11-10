# New method implementations for generic functions
# 
# Author: Peter Humburg
###############################################################################

####################### methods for strandPileup ##############################
## strand specific read counts from 'AlignedRead'
setMethod("strandPileup", "AlignedRead",
		definition=function(aligned, chrLen, extend, coords=c("leftmost", "fiveprime"), 
				compress=TRUE, plot=TRUE, ask=FALSE, ...){
			
			coords <- match.arg(coords)
			chrFilter <- lapply(levels(chromosome(aligned)), function(chr, all) all == chr, chromosome(aligned))
			
			## compute pileup for each chromosome
			counts <- mapply(function(filter, len, extend, aligned, ...){ 
#						sfilter <- filter & strand(aligned) == "+"
#						start1 <- position(aligned)[sfilter]
#						sfilter <- filter & strand(aligned) == "-"
#						start2 <- position(aligned)[sfilter]
#						width2 <- width(aligned)[sfilter]
#						counts <- cbind(pileup(start1, fraglength=extend, chrlength=len, 
#										factor("+", levels=c("-", "+", "*")), ...),
#								pileup(start2, fraglength=extend, chrlength=len, 
#										factor("-", levels=c("-", "+", "*")), 
#										readlength=width2-1, ...))
						sfilter1 <- filter & strand(aligned) == "+"
						sfilter2 <- filter & strand(aligned) == "-"
						## separate strands
						fwd <- aligned[sfilter1]
						rev <- aligned[sfilter2] 
						ext1 <- -width(fwd) + as.integer(extend)
						ext2 <- -width(rev) + as.integer(extend)
						counts <- list(coverage(fwd, extend = ext1, coords = coords, ...)[[1]],
								coverage(rev, extend = ext2, coords = coords, ...)[[1]])
						## append trailing zeros to get full chromosome length
						counts <- .fixCounts(counts, len)
						
						if(!compress) counts <- decompress(counts)
						counts
					},  
					chrFilter, chrLen, extend, MoreArgs=c(list(aligned=aligned), list(...)), SIMPLIFY=FALSE)
			
			counts <- ReadCounts(counts, levels(chromosome(aligned)), compress=compress)
			if(!plot) return(counts)
			
			## plot image of counts for each chromosome and strand 
			ask <- devAskNewPage(ask)
			for(i in 1:length(counts)){
				plot(counts, chr=i, type="hilbert", scale="ratio", log=TRUE)
				plot(counts, chr=i, type="hilbert", scale="total", log=TRUE)
			}
			devAskNewPage(ask)
			
			invisible(counts)
		}
)

## strand specific read counts from data frame
setMethod("strandPileup", "data.frame",
		definition=function(aligned, chrLen, extend, coords=c("leftmost", "fiveprime"), 
				compress=TRUE, plot=TRUE, ask=FALSE, ...){
			
			coords <- match.arg(coords)
			## ensure aligned contains all required information
			## Note that all operations are based on column names. This allows
			## for the presence of additional columns and arbitrary column order
			## but imposes strict requirements on column names.
			requiredNames <- c("chromosome", "start", "end", "strand")
			
			## need names for chrLen
			if(is.null(names(chrLen))) names(chrLen) <- levels(aligned$chromosome)
			
			## if we have start position and length for each read we convert this
			## into start and end position
			if(all(c("position", "length") %in% names(aligned)) && !isTRUE("start" %in% names(aligned)))
				names(aligned)[which(names(aligned) == "position")] <- "start"
			## ensure 'start' is 5'-end and 'end' is 3'-end
			if(all(c("start", "length") %in% names(aligned)) && !isTRUE("end" %in% names(aligned))){
				idx <- aligned$strand == "+"
				if(coords == "leftmost") {
					aligned$end <- aligned$start + aligned$length - 1
					tmp <- aligned$start[!idx]
					aligned$start[!idx] <- aligned$end[!idx]
					aligned$end[!idx] <- tmp
				}
				if(coords == "fiveprime"){
					aligned$end <- integer(nrow(aligned))
					aligned$end[idx] <- aligned$start[idx] + aligned$length[idx] - 1
					aligned$end[!idx] <- aligned$start[!idx] - aligned$length[!idx] + 1
				}		
			}
			if(!all(requiredNames %in% names(aligned)))
				stop("Columns ", requiredNames[which(!(requiredNames %in% names(aligned)))], " are missing.")
			
			## compute pileup for each chromosome
			counts <- vector(length(levels(aligned$chromosome)), mode="list")
			names(counts) <- levels(aligned$chromosome)
			for(chr in levels(aligned$chromosome)){
				chrFilter <- aligned$chromosome == chr
				strFilter <- aligned$strand == "+"
				start1 <- aligned$start[chrFilter & strFilter]
				start1 <- start1[start1 + extend -1 <= chrLen[chr]]
				start2 <- aligned$end[chrFilter & !strFilter] - extend + 1
				start2 <- start2[start2 >= 1L]
				
#				if(extend > 1){
#					start1 <- unlist(lapply(start1, function(x) seq(x, x+extend-1)))
#					start2 <- unlist(lapply(start2, function(x) seq(x, x-extend+1)))
#				} 
#				counts[[chr]] <- cbind(pileup(aligned$start[chrFilter & strFilter], fraglength=extend, 
#								chrlength=chrLen, dir=factor("+",levels=c("-", "+", "*")), ...),
#						pileup(aligned$end[chrFilter & !strFilter], fraglength=extend, chrlength=chrLen, 
#								readlength=extend-1,dir=factor("-",levels=c("-", "+", "*")), ...))
				
				counts[[chr]] <- list(coverage(IRanges(start=start1, width=extend), ...), 
						coverage(IRanges(start=start2, width=extend), ...))
				counts[[chr]] <- .fixCounts(counts[[chr]], chrLen[chr])
				if(!compress) counts[[chr]] <- decompress(counts[[chr]])
			}
			counts <- ReadCounts(counts, compress=compress)
			if(!plot) return(counts)
			
			## plot image of counts for each chromosome and strand
			ask <- devAskNewPage(ask)
			for(i in 1:length(counts)){
				plot(counts, chr=i, type="hilbert", scale="ratio", log=TRUE)
				plot(counts, chr=i, type="hilbert", scale="total", log=TRUE)
			}
			devAskNewPage(ask)
			
			invisible(counts)
		}
)


####################### methods for callBindingSites ##############################

## binding sites from list of read counts (one component per chromosome) 
## this is where the actual work is done, other methods should create strand specific 
## read counts and pass them on to this method
setMethod("callBindingSites", "ReadCounts", 
		definition=function(data, bind, support, background, bgCutoff=0.9, supCutoff=0.9, 
				fdr = 0.05, extend=1, tailCut=0.95, piLambda=0.5, adapt=FALSE, corSummary=median, 
				compress=TRUE, digits=16, plot=TRUE, verbose=TRUE, ask=FALSE, plotTo, ...){
			
			if(missing(bind))
				stop("Mandatory argument \"bind\" is missing. Please specify length of binding region. ",
						"You may specify a range of possible values by providing the minimum and maximum.")
			if(missing(support))
				stop("Mandatory argument \"support\" is missing. Please specify length of support region. ",
						"You may specify a range of possible values by providing the minimum and maximum.")
			
			
			## plot to file plotTo
			if(!missing(plotTo) && plot) pdf(plotTo, width=8, height=11)
						
			## determine length of binding and support region if necessary
			if(length(bind) > 1 || length(support) > 1){
				if(verbose) message("Determining length of binding and support region...")
				if(plot) par(mfrow=c(2,1))
				regionLen <- getBindLen(data, bind, support, corSummary, verbose=verbose, plot=plot)
				bind <- regionLen[1]
				support <- regionLen[2]
			}
			
			## ensure 'background' has sensible value
			if(missing(background)){
				background <- 10*(bind + 2*support)
				if(verbose) message("Setting length of background window to ", background, ".")
			}
			if(background < 2*(bind + 2*support)){
				warning("Background window is too short (", background, "). Using", 10*(bind + 2*support), "instead.\n")
				background <- 10*(bind + 2*support)
			}
			
			## calculate score
			if(verbose) message("Scoring binding sites...")
			score <- lapply(data, startScore, b=bind, support=support, background=background, bgCutoff=bgCutoff, 
					supCutoff=supCutoff)
			
			## identify significant peaks
			if(verbose) message("Determining significance...")
			allScores <- c(lapply(score, function(x) x[!is.na(x)]), recursive=TRUE)
			cutoff <- getCutoff(allScores, alpha = fdr, lambda = piLambda, tailCut=tailCut, adapt=adapt, plot = plot, 
					returnPval=TRUE)
			if(verbose) message("Significance threshold [score (p-value)]: ", sprintf("%.2f", cutoff$cutoff[1]), " (", 
						cutoff$cutoff[2], ")")
			rm(allScores)
			peaks <- lapply(score, pickPeak, cutoff$cutoff[1], offset=support + ceiling(bind/2))
			
			
			## need to map p-values back to genomic coordinates
			at <- lapply(score, is.na)
			pval <- vector(length(score), mode="list")
			start <- 0
			for(i in 1:length(score)){
				pval[[i]] <- numeric(length(score)) + NA
				pval[[i]][which(!at[[i]])] <- cutoff$pvalue[1:sum(!at[[i]])+start]
				start <- start + sum(!at[[i]])
			}
			
			call <- match.call(callBindingSites, call=sys.call(sys.parent()))
			call$bind <- eval(call$bind)
			call$support <- eval(call$support)
			call$background <- eval(call$background)
			start <- support + ceiling(bind/2)
			result <- BindScore(call=call, score=score, pvalue=pval, peaks=peaks, cutoff=cutoff$cutoff, 
					nullDist=cutoff$h0, names=names(data), compress=compress, digits=digits, start=start)
			
			if(plot){
				par(mfrow=c(1, 2))
				plot(result, type="density")
				plot(result, type="qqplot")
			}
			if(!missing(plotTo) && plot) dev.off()
			
			result
		}
)

## using matrix of read counts (single chromosome)
setMethod("callBindingSites", "matrix",
		definition=function(data, chrName="chr", ...){
			data <- new("ReadCounts", list(data), chrName)
			call <- match.call(callBindingSites,sys.call())
			call$data <- data
			call[[1]] <- `callBindingSites`
			for(i in 3:length(call)) call[[i]] <- eval(call[[i]])
			eval(call)
		}
)

## using name of file with mapped reads
## type: file type (see ?readAligned for details)
## minQual: minimum alignment quality to use. All reads with lower alignment quality will be removed 
##          before pilup is created
setMethod("callBindingSites", "character",
		definition=function(data, type, minQual=70, ...){
			data <- readAligned(data, type=type)
			data <- data[quality(alignQuality(data)) >= minQual]
			
			call <- match.call(callBindingSites,sys.call())
			call$data <- data
			call[[1]] <- `callBindingSites`
			for(i in 3:length(call)) call[[i]] <- eval(call[[i]])
			eval(call)			
		})


## For everything else we try to compute read counts and pass data on to the list method
setMethod("callBindingSites", "ANY",
		definition=function(data, chrLen, plot=TRUE, verbose=TRUE, ..., plotTo){
			## plot to file plotTo
			if(!missing(plotTo) && plot) {
				pdf(plotTo, width=8, height=11)
				par(mfrow=c(length(chrLen), 2))
			}
			if(verbose) message("Accumulating read counts...")
			
			dots <- list(...)
			if(is.null(dots$extend)) dots$extend <- 1			 
			counts <- strandPileup(data, chrLen=chrLen, extend=dots$extend, plot=plot)
			
			call <- match.call(callBindingSites,sys.call())
			call$data <- counts
			call[[1]] <- `callBindingSites`
			for(i in 3:length(call)) call[[i]] <- eval(call[[i]])
			bindScore <- eval(call)
			
			if(!missing(plotTo) && plot) dev.off()
			
			return(bindScore)
		}
)

############################# New methods for existing generic functions ##############################

## lapply method for ReadCounts objects
setMethod("lapply", "ReadCounts",
		definition=function(X, FUN, ...){
			lapply(X@counts, FUN, ...)
		}
)

## sapply method for ReadCounts objects
setMethod("sapply", "ReadCounts",
		definition=function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE){
			sapply(X@counts, FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES)
		}
)

## lapply for BindScore objects. Applies function by chromosome
setMethod("lapply", "BindScore",
		definition=function(X, FUN, ...){
			result <- vector(length(names(X)), mode="list")
			names(result) <- names(X)
			for(i in names(X)){
				result[[i]] <- match.fun(FUN)(X[[i]], ...)
			}
			result
		}
)


## number of chromosomes
setMethod("length", "ReadCounts", definition=function(x) length(x@counts))

## length of chromosomes
setMethod("chrLength", "ReadCounts", definition=function(x, subset){
			len <- sapply(x, nrow)
			if(!missing(subset)) len <- len[subset]
			len
		} 
)

setMethod("chrLength", "RLEReadCounts", definition=function(x, subset){
			len <- sapply(x, IRanges::elementLengths)[1, ]
			if(!missing(subset)) len <- len[subset]
			len
		} 
)

setMethod("chrLength", "BindScore", definition=function(x, subset){
			len <- IRanges::elementLengths(score(x)) + 2*support(x) + binding(x)
			if(!missing(subset)) len <- len[subset]
			len
		} 
)

## number of reads on each chromosome/strand
setMethod("nreads", "ReadCounts", definition=function(x, byStrand = TRUE, subset){
			if(byStrand){ 
				counts <- sapply(x, colSums)
				if(!missing(subset)) counts <- matrix(counts[ , subset], nrow = 2, 
							dimnames=list(c("+", "-"), colnames(counts)[subset]))
			} else{ 
				counts <- sapply(x, sum)
				if(!missing(subset)) counts <- counts[subset]
			}
			counts
		}
)

setMethod("nreads", "RLEReadCounts", definition=function(x, byStrand = TRUE, subset){
			counts <- sapply(x, sum)
			if(!missing(subset))
				counts <- matrix(counts[ , subset], nrow = 2, dimnames=list(c("+", "-"), 
								colnames(counts)[subset]))
			if(!byStrand) counts <- colSums(counts)
			counts
		}
)


## number of predicted binding sites
setMethod("length", "BindScore", definition=function(x) sum(sapply(peaks(x), length)))

setMethod("length<-", "ReadCounts",
		definition=function(x, value){length(x@counts) <- value; invisible(x)})

## TODO: this is inconsistent with the length method above, consider changing
setMethod("length<-", "BindScore", 
		definition=function(x, value){
			length(x@score) <- value
			length(x@pvalue) <- value
			length(x@peaks) <- value
			invisible(x)
		}
)

## convert to data.frame
setMethod("as.data.frame", c(x="BindScore"),
		definition=function(x, ...){
			chr <- rep(names(x), times=IRanges::elementLengths(peaks(x)))
			pos <- c(peaks(x), recursive=TRUE)
			peakScores <- c(mapply(function(y, i, offset) as.numeric(y[i - offset]),  score(x), peaks(x), 
							MoreArgs=list(support(x) + ceiling(binding(x)/2)), SIMPLIFY=FALSE), 
					recursive=TRUE)
			peakPval <- c(mapply(function(y, i, offset) as.numeric(y[i - offset]), pvalue(x), peaks(x), 
							MoreArgs=list(support(x) + ceiling(binding(x)/2)), SIMPLIFY=FALSE), 
					recursive=TRUE)
			
			df <- data.frame(chromosome=chr, position=pos, score=peakScores, pvalue=peakPval)
			if(nrow(df) > 0) rownames(df) <- 1:nrow(df)
			
			df
		}
)


## get first n peaks (either by position or score)
setMethod("head", "BindScore", 
		definition=function(x, n=6, by=c("score", "position"), ...){

			df <- as.data.frame(x)
			
			type <- pmatch(by[1], c("score", "position"), nomatch=0)
			if(type == 1){
				idx <- order(df$pvalue, -df$score)
				df <- df[idx,]
			} 
			
			head(df, n, ...)
		}
)

## get last n peaks (either by position or score)
setMethod("tail", "BindScore", 
		definition=function(x, n=6, by=c("score", "position"), ...){
			
			df <- as.data.frame(x)
			
			type <- pmatch(by[1], c("score", "position"), nomatch=0)
			if(type == 1){
				idx <- order(df$pvalue, -df$score)
				df <- df[idx,]
			} 
			
			tail(df, n, ...)
		}
)

## minimum and maximum score
setMethod("min", "BindScore",
		definition=function(x, ..., na.rm=TRUE){
			min(sapply(score(x), min, na.rm=na.rm))
		}
)

setMethod("max", "BindScore",
		definition=function(x, ..., na.rm=TRUE){
			max(sapply(score(x), max, na.rm=na.rm))
		}
)

setMethod("range", "BindScore",
		definition=function(x, ..., na.rm=TRUE){
			c(min(x, na.rm=na.rm), max(x, na.rm=na.rm))
		}
)

## human readable summaries
setMethod("show", "ReadCounts", 
		definition=function(object){
			.showHeader(object)
			cat("Number of reads:", sum(sapply(object@counts, sum)),"\n")			
		}
)

setMethod("show", "BindScore", 
		definition=function(object){
			.showHeader(object)
			cat("Predicted binding sites:", length(object), "\n")
			cat("Significance cut-off: ", cutoff(object, "score"), 
					" (p-value: ", cutoff(object, "pvalue"), ")\n", sep = "")
		}
)

## explicit conversion between compressed and decompressed representations
setMethod("compress", "ReadCounts",
		definition=function(x){
			ReadCounts(counts=x@counts, compress=TRUE)
		}
)
setMethod("compress", "BindScore",
		definition=function(x, digits=16){
			BindScore(call=x@functionCall, score=score(x), pvalue=pvalue(x), peaks=peaks(x),
					cutoff=cutoff(x), nullDist=nullDist(x), compress=TRUE, digits=digits)
		}
)
setMethod("compress", "RLEReadCounts",
		definition=function(x){ ## input is already compressed, do nothing
			x
		}
)
setMethod("compress", "RLEBindScore",
		definition=function(x){ ## input is already compressed, do nothing
			x
		}
)
setMethod("decompress", "RLEReadCounts", 
		definition=function(x){
			ReadCounts(counts=lapply(x@counts, decompress, simplify=TRUE), compress=FALSE)
		}
)
setMethod("decompress", "RLEBindScore", 
		definition=function(x){
			BindScore(call=x@functionCall, score=lapply(score(x), as.numeric), 
					pvalue=lapply(pvalue(x), as.numeric), peaks=peaks(x), cutoff=cutoff(x), 
					nullDist=nullDist(x),start=x@start, compress=FALSE)
		}
)

setMethod("decompress", "Rle",
		definition=function(x, class){
			if(missing(class)) class <- base::class(runValue(x))
			as(x, class)
		}
)
setMethod("decompress", "RleList",
		definition=function(x, class, simplify=TRUE){
			if(missing(class)) class <- sapply(x, function(i) base::class(runValue(i)))
			result <- lapply(1:length(x), function(i) as(x[[i]], class[i]))
			names(result) <- names(x)
			
			## simplify into matrix or vector if possible
			if(simplify & length(result) == 1) result <- result[[1]]
			else if(simplify & all.equal(class, rep(class[1], length(class)), check.attributes=FALSE) & 
					max(abs(diff(IRanges::elementLengths(result)))) == 0){
				result <- do.call(cbind, result)
			}
			
			result
		}
)

#################  plotting functions ##########################

## plot read counts in a window, optionally with overlayed score and nucleosome positions
## or plot hilbert curve of all read counts on chromosome
setMethod("plot", list(x="ReadCounts", y="missing"), 
		definition=function(x, chr, center, score, width=2000, type=c("hilbert", "window"), ...){
			type <- match.arg(type)
			if(missing(chr) && length(x) == 1) chr <- names(x)
			if(type == "window"){
				if(!missing(score)){
					start <- center - ceiling(width/2)
					score <- decompress(score[[chr, start:(start + width+1)]])
				} 
				.plotWindow(x, chr=chr, center=center, score=score, width=width, ...)
			}
			else if(type == "hilbert") .plotReads(x[[chr]], ...)
		}
)

setMethod("plot", list(x="RLEReadCounts", y="missing"), 
		definition=function(x, chr, center, score,  width=2000, type=c("hilbert", "window"), ...){
			type <- match.arg(type)
			if(missing(chr) && length(x) == 1) chr <- names(x)
			if(type == "window"){
				start <- floor(center - ceiling(width/2))
				end <- ceiling(center + ceiling(width/2))
				if(!is.character(chr)) chr <- names(x)[chr]
				x <- list(sapply(x[[chr]], function(x) as.integer(window(x,start,end))))
				names(x) <- chr
				
				if(!missing(score)) score <- decompress(score[[chr, start:end]])
				.plotWindow(x, chr=chr, start=start, end=end, width=width, score=score, offset=start, ...)
			}
			else if(type == "hilbert"){
				.plotReads(decompress(x[[chr]]), ...)
			}
		}
)

## plots to asses fit of null distribution
setMethod("plot", list(x="BindScore", y="missing"),
		definition=function(x, npoints = 10000, type=c("density", "qqplot"), ...){
			type <- match.arg(type)
			
			cutoff <- cutoff(x)["score"]
			nullDist <- nullDist(x)
			## get vector of all scores
			values <- unlist(lapply(score(x), function(x) {d <- decompress(x); d[!is.na(d)]}),
					use.names=FALSE)
			if(type == "density"){
				minVal <- min(values)
				maxVal <- max(values)
				range <- maxVal - minVal
				breaks <- seq(minVal, maxVal, length.out=max(range, 100))
				hist(values, breaks=breaks, freq=FALSE, main="Empirical and fitted null distribution", xlab="Score")
				lines(seq(minVal, maxVal, by=0.1), dnorm(seq(minVal, maxVal, by=0.1),
								mean=nullDist["mean"], sd=nullDist["sd"]), col=2)
				if(is.finite(cutoff)) abline(v=cutoff, col=4, lty=2)
			}
			if(type == "qqplot"){
				## qqnorm is very slow for large samples
				## if the sample is large we just plot a subset of all points
				if(length(values) <= npoints){ 
					qqnorm(values)
					qqline(values, col=2)
				}
				else{
					qObs <- quantile(values, probs=seq(0,1, length.out=npoints))
					qTheo <- qnorm(ppoints(qObs))
					plot(qTheo, qObs, xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="Normal Q-Q Plot")
					qqline(values, col=2)
				}
				if(is.finite(cutoff)) abline(h=cutoff, col=4, lty=2)
			}
		})