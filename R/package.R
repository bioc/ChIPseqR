# Functions for ChIPseqR package
# 
# Author: Peter Humburg
###############################################################################


############## read mapping #######################

## align sequence data according to transcription start sites
## data: mapped reads (ReadCounts object)
## anno: gff formated annotations for genes (or any other feature of interest)
## offset: number of base pairs to use from either side of TSS
## returns a list of sequences centred around TSS
alignFeature <- function(data, anno, offset = 1000){
	align <- function(start,x,offset){
		
		chr <- as.character(start[1,1])
		pos <- as.numeric(start[1,2])
		
		if(is.null(x[[chr]])) return(NULL)
		
		## read counts may be compressed
		x <- decompress(x[[chr]])
		
		## may be dealing with pooled forward and reverse read counts
		if(!is.matrix(x)) x <- matrix(x, ncol=1)
		
		ret <- matrix(0,ncol=ncol(x),nrow=2*offset+1)
		slice <- (offset-min(pos,offset)+1):(offset+min(offset,nrow(x)-pos-1)+1)
		ret[slice,] <- x[(pos-(offset-slice[1]+1)):(pos+(slice[length(slice)]-offset-1)),]
		
		ret
	}
	seq.genes <- vector(nrow(anno), mode="list")
	
	for(i in 1:nrow(anno)){
		reverse <- anno[i,7]=="-"
		seq.genes[[i]] <- align(anno[i,c(1,4+reverse)], data, offset)
		if(reverse) seq.genes[[i]] <- rev(seq.genes[[i]])
	}
	
	
	seq.genes	
}

## compute sliding window summaries of read counts
windowCounts <- function(reads, window=1000, shift=500, method=sum){
	n <- length(reads)
	start <- seq(1, n-window, by=shift)
	end <- seq(window, n, by=shift)
	wndw <- aggregate(reads, start=start, end=end, FUN=method)
	names(wndw) <- ceiling(start+(end-start)/2)
	wndw
}


########### visualisation ####################


## plot read counts in a window, optionally with overlayed score and nucleosome positions
.plotWindow <- function(data, chr, center, score, width=2000, bind, start, end,
		bind.col=3, score.type='l',	xlab=NULL, ylab="Read count", cutoff=TRUE, offset=1, ...){
	stopifnot(is(data, "ReadCounts") || is.list(data))
	
	data <- data[[chr]]
	if(missing(start)) start <- floor(center - ceiling(width/2))
	if(missing(end)) end <- ceiling(center + ceiling(width/2))
	
	range <- start:end - offset + 1
	height <- max(data[range,])
	
	if(!missing(score)){
		mar <- par("mar")
		mar[4] <- max(4.1, mar[4])
		par(mar=mar)
		score <- score[[chr, start:end]]
		peaks <- unlist(peaks(score))
		if(missing(bind)) bind <- binding(score)
		support <- support(score)
		binding <- binding(score)
		cutValue <- cutoff(score)["score"]
		score <- unlist(score(score))
	}
	
	## plot read counts
	if(is.null(xlab)) xlab <- paste(chr, start, "-", end)
	plot(start:end, data[range,1], ylim=c(0,height), type='h', ylab=ylab, xlab=xlab, col=2, ...)
	lines(start:end, data[range,2], type='h', col=4, lty=2)
	
	## plot score
	if(!missing(score)){
		score.col <- col2rgb(3)/255
		score.col <- rgb(sqrt(score.col[1,]), sqrt(score.col[2,]), sqrt(score.col[3,]))
		
		min.score <- min(score, na.rm=TRUE)
		max.score <- max(score, na.rm=TRUE)
		denom <- ((max.score-min.score)/height)
		lines(start:end, (score - min.score)/denom, col=score.col, lty=1, type=score.type, cex=0.5)
		if(cutoff) abline(h=(cutValue - min.score)/denom, col=score.col, lty=2)
		
		## add axis
		ticks <- pretty(score)
		ticks.at <- height*(ticks - min.score)/(max.score-min.score)
		axis(4, labels=ticks, at=ticks.at, col=score.col, col.axis=score.col)
		mtext("Score", side=4, line=2.5, cex=par("cex")*par("cex.lab"), col=score.col)
		
		## plot predicted binding sites
		if(length(peaks) > 0)
			.plotBind(peaks, 0, col=bind.col, b=bind, bind=binding, offset=0, 
					step=-0.01*(par("usr")[4]-par("usr")[3]), lwd=4, lend=1, extend=support)
	}
	
}

## plot binding sites
## plots binding region (bind) of predicted binding site as well as actual (assumed) length of binding site (b) 
.plotBind <- function(x, y, col=1, add=TRUE, b=147, bind=128, offset=-1000, step, extend, ...){
	## assume bind <= b
	## TODO: May want to allow bind > b, not sure whether it is needed in practice though
	x.centre <- x
		
	## calculate vertical position of features to plot
	len <- max(b, bind + 2*extend) + 1
	ypos <- numeric(length(x)) + y
	if(length(x) > 1) for(i in 2:length(x)){
		if(x[i] <= x[i-1] + len) ypos[i] <- ypos[i-1] + step
	}
	
	x <- x.centre - ceiling(b/2)
	## plot larger feature first
	if(!identical(b, bind)){
		col <- rep(col,length.out=length(x))
		if(!add){
			plot(c(x[1],x[1]+b-1)+offset, c(y,y), col=col[1], type='l',...)
		}
		
		for(i in (1+!add):length(x)){
			null<-mapply(function(s,col){ lines(c(s,s+b)+offset,c(ypos[i],ypos[i]),col=col,...)}, x[[i]],col[[i]])
		}
		add <- TRUE
	}
	
	## use lighter colour for overlayed feature
	col <- col2rgb(col)/255
	col <- rgb(sqrt(col[1,]), sqrt(col[2,]*0.5), sqrt(col[3,]))
	
	x <- x.centre - ceiling(bind/2)
	col <- rep(col,length.out=length(x))
	if(!add){
		plot(c(x[1],x[1]+bind-1)+offset, c(y,y), col=col[1], type='l',...)
		if(!missing(extend) && extend > 0){
			lines(c(x-extend, x)+offset, c(y,y), col=col, lend=1)
			lines(c(x + bind + extend - 1, x)+offset, c(y,y), col=col, lend=1)
		} 
	}
	
#	pos <- !add
	for(i in (1+!add):length(x)){
#		if(i > 1 && x[[i]] > x[[i-1]] + bind + ifelse(missing(extend), 0, 2*extend)) pos <- 0
		null<-mapply(function(s,col){ lines(c(s,s+bind)+offset,c(ypos[i],ypos[i]),col=col,...)}, x[[i]],col[[i]])
		if(!missing(extend) && extend > 0){
			null <- mapply(function(x,col) lines(c(x-extend, x)+offset, c(ypos[i],ypos[i]), col=col, lend=1), 
					x[[i]], col[[i]])
			null <- mapply(function(x,col) lines(c(x + bind + extend, x)+offset, c(ypos[i],ypos[i]), col=col, lend=1), 
					x[[i]], col[[i]])
		}
#		pos <- pos + 1
	}
}

.plotReads <- function(x, scale=c("total", "ratio"), log=TRUE, ...){
	stopifnot(is(x, "matrix"), ncol(x) == 2)
	method <- match.arg(scale)
		
	hilbert <- list(HilbertVis::hilbertImage(x[,1]), HilbertVis::hilbertImage(x[,2]))
	
	## process plot parameters
	dots <- list(...)
	if(!is.null(dots$mar)) mar <- dots$mar
	else mar <- c(0,0,4,0)
	if(!is.null(dots$col)) col <- dots$col
	
	oldPar <- par("mar","bg","xpd")
	
	par(mar=mar, bg="gray", xpd=TRUE)
	if(method == "total"){
		data <- hilbert[[1]] + hilbert[[2]]
		if(log) data <- log(data)
		
		if(is.null(dots$col)) col <- fBasics::seqPalette(256,"Blues")
		dataRange <- diff(range(data[is.finite(data)], na.rm=TRUE))
		centre <- as.numeric(!log)
		step <- dataRange/length(col)
		breaks <- c(seq(min(data[is.finite(data)], na.rm=TRUE), centre - step, length.out=length(col)/2), centre, 
				seq(step, max(data[is.finite(data)], na.rm=TRUE), length.out=length(col)/2))
		
		image(data, col=col, breaks=breaks, xaxt="n",yaxt="n")
	}
	else if(method == "ratio"){
		if(log)
			data <- log(hilbert[[1]]) - log(hilbert[[2]])
		else data <- hilbert[[1]]/hilbert[[2]]
		if(is.null(dots$col)) col <- fBasics::divPalette(256,"RdYlGn")
		image(data, col=col,xaxt="n",yaxt="n")		
	}
	## ensure margin background is white
	offset <- 1/(2*(nrow(data)-1)) 
	top <-  grconvertY(grconvertY(1 + offset, "user", "inches") + par("mai")[3], "inches", "user")
	rect(-offset, 1 + offset, 1 + offset, top, col="white", border=NA)
	
	## add title
	if(!is.null(dots$main)) title(dots$main)
	else{
		if(method == "ratio") main <- "Ratio of Read Counts"
		else if(method == "total") main <- "Total Read Counts"
		title(paste(ifelse(log, "Log", ""), main))	
	}
	
	par(oldPar)	
}

####################### conversion #####################################
## convert list of genome coordinates into gff format
## pos: named list with start positions for each chromosome
pos2gff <- function(pos, method, feature, len, strand, score, name){
	total <- sum(IRanges::elementLengths(pos))
	chr <- rep(names(pos), times = IRanges::elementLengths(pos))
	start <- unlist(pos, use.names=FALSE)
		
	method <- rep(method, length.out=total)
	feature <- rep(feature, length.out=total)
	if(missing(strand)) strand <- "."
	strand <- rep(strand, length.out=total)
	if(missing(score)) score <- "."
	score <- rep(score, length.out=total)
	if(missing(name)) name <- paste(feature, 1:total, sep="_")
	name <- paste("name \"", name, "\"", sep='')
	
	data.frame(chromosome=chr, method=method, feature=feature, start=start, 
			end=start+len, score=score, strand=strand, frame=".", name=I(name),
			stringsAsFactors=TRUE)
}


####################### finding binding sites ##########################

## using C code
startScore <- function(data, b, support, background, bgCutoff, supCutoff){
	if(typeof(data) != "integer")
            data <- matrix(as.integer(unlist(data, use.names=FALSE)), ncol = 2)
	score <- .Call(startScore_pois, data, as.integer(b), as.integer(support), as.integer(background), 
			bgCutoff, supCutoff)
	score
}

## given binding site scores (for one chromosome), identify all peaks above given threshold
pickPeak <- function(score, threshold, offset=0, sub=FALSE){
	idx <- which(score > threshold)
	if(length(idx) == 0) return(NULL)
	
	idx.dist <- diff(idx)
	change.at <- which(idx.dist > 1)
	idx.start <- c(1, change.at +1)
	idx.end <- c(change.at, length(idx))
	broadPeak <- IRanges::Views(score, idx[idx.start], idx[idx.end])
	peak <- IRanges::viewWhichMaxs(broadPeak)
	## if the peak is flat we pick the centre
	for(i in 1:length(peak)){
	  runs <- as(broadPeak[[i]][(peak[i] - start(broadPeak)[i] + 1):width(broadPeak)[i]], "Rle")
	  peakWidth <- runLength(runs)[1]
	  if(peakWidth > 1){
	    peak[i] <- peak[i] + floor(peakWidth/2) 
	  }
	}
	if(sub){
		subPeaks <- mapply(function(s,e, x) {
					if(s == e) p <- s + offset
					else if(isTRUE(all.equal(diff(IRanges::window(x, s, e)), rep(0, e-s)))) 
						p <- s + floor((e - s)/2) + offset
					else if(e-s > 1){
						p <- which(diff(IRanges::window(x, s, e-1)) >= 0 &  
										diff(IRanges::window(x, s+1, e)) <= 0) + s + offset
						if(x[s] > x[s + 1]) p <- c(s + offset, p)
						if(x[e] > x[e - 1]) p <- c(p, e + offset)
						drop <- which(diff(p) == 1)
						if(length(drop) > 0) p <- p[-(drop + 1)]
					}
					else p <- which.max(x[s:e]) + s -1 + offset
					p
				},idx[idx.start], idx[idx.end], 
				MoreArgs=list(x=score), SIMPLIFY=FALSE)
		ret <- list(peaks=peak+offset, subPeaks=subPeaks)
	}
	else ret <- peak + offset
	
	return(ret)
}

## determine significance threshold for binding site scores
## alpha is the target level of the FDR (which will be adapted to the data)
getCutoff <- function(score, alpha = 0.05, tailCut=0.95, adapt=FALSE, lambda, plot = TRUE, returnPval = TRUE){
	nllk <- function(sd, x, cutoff){
		dnom <- log(pnorm(cutoff, sd=sd, log.p=FALSE) - pnorm(0, sd=sd, log.p=FALSE))
		-(sum(dnorm(x[x < cutoff], sd=sd, log=TRUE) - dnom))
	}
	
	score <- score[!is.na(score)]
	score2 <- -score[score <= median(score)]
	score2 <- score2 - min(score2)
	
	fit <- optimize(nllk, c(.Machine$double.eps, .Machine$double.max), x=score2, cutoff=quantile(score2, tailCut))
	sigma <- fit$minimum
	pval <- pnorm(score, mean=median(score), sd=sigma, lower.tail=FALSE)
	fdr <- p.adjust(pval, "BH")
	
	pi0 <- NA
	if(adapt){
		F_hat <- sum(pval <= lambda)/length(pval)
		pi0 <- min(1, (1 - F_hat)/(1 - lambda))
		alpha <- alpha/pi0
	}
	
	cutoff <- ifelse(min(fdr) > alpha, Inf, min(score[fdr <= alpha])) 
	
#	if(plot){
#		range <- max(score) - min(score)
#		breaks <- seq(min(score), max(score), length.out=max(range, 100))
#		hist(score, breaks=breaks, freq=FALSE, main="Empirical and fitted null distribution", xlab="Score")
#		lines(seq(min(score), max(score), by=0.1), dnorm(seq(min(score), max(score), by=0.1),
#						mean=median(score), sd=sigma), col=2)
#		if(is.finite(cutoff)) abline(v=cutoff, col=4, lty=2)
#	}
	
	## assemble result
	result <- list(cutoff=c(cutoff, alpha), h0=c(median(score), sigma))
	names(result$cutoff) <- c("cutoff", "alpha")
	names(result$h0) <- c("mean", "sd")
	
	if(returnPval) result[["pvalue"]] <- fdr
	if(adapt) result[["pi0"]] <- pi0
	
	return(result)
}

## calculate cross-correlation between strands
getBindCor <- function(data, max.lag, summary, plot=TRUE, ...){
	result <- sapply(lapply(decompress(data), timsac::fftcor, lag=max.lag, plot=FALSE), function(x) x$ccor21)
	
	if(plot){
		matplot(0:max.lag,result, type='l', xlab="Lag", ylab="Cross-correlation", 
				main="Correlation between forward and reverse strand", lty=1:5, col=1:6)
		if(is.null(names(data))) names <- 1:length(data)
		else names <- names(data)
		if(length(names) > 1 )legend("topright", legend=names, col=1:6, lty=1:5)
	}
	
	if(!missing(summary)) result <- apply(result, 1, summary)
	
	result
}

## identify peaks in cross-correlation between strands to determine binding site length (and optionally support region length)
## bind and support give the maximum (and optionally minimum) length for binding site and support region
##
## Note that the assumption that the first peak in the cross-correlation indicates the length of the binding site 
## is not accurate. The peak is closer to bind + 2*m where m is the median of the read distribution in the support region.
## ('read distribution in the support region' means the read density as a function of distance to binding site start/end)
## Consequently this method will overestimate the length of the binding site.
## If either bind or support are of length 1 this is assumed to be the known value and a more accurate estimate for the
## remaining parameter is used.
getBindLen <- function(data, bind, support, summary=median, verbose=FALSE, plot=TRUE, ...){
	max.lag <- ifelse(missing(support), ceiling(max(bind)*1.5), 2*max(bind)+ceiling(max(support)*1.5))
	bindCor <- getBindCor(data, max.lag, summary, plot=FALSE, ...)
	
	## fit smoothing spline 
	bindSpline <- smooth.spline(0:(length(bindCor)-1), bindCor)
	
	#determine lag that maximises correlation, this should be close to length of binding site
	if(length(bind) == 1){
		bindLen <- bind
		#supLen <- which.max(bindSpline$y[(bind + min(support)):(bind + max(support)) + 1]) + min(support) - 1
		supPeak <- which.max(IRanges::window(bindSpline$y, bindLen + 2*min(support),
						bindLen + 2*max(support) + 1)) + bindLen + 2*min(support) - 1
		supLen <- round((supPeak - bindLen) * 0.5)
		bindPeak <- supPeak   ## need only location of first peak 
		if(verbose) message("Estimated length of support region: ", supLen)
	}
	else{
		if(length(support) == 1){
			supLen <- support
			bindPeak <- which.max(IRanges::window(bindSpline$y, supLen + min(bind), 
									supLen + max(bind) + 1)) +	min(bind) + supLen -1 
			bindLen <- bindPeak - 2*supLen
			supPeak <- bindPeak   ## need only location of first peak
			if(verbose) message("Estimated length of binding site: ", bindLen)
		}
		else{
			bindPeak <- which.max(IRanges::window(bindSpline$y, support[1] + bind[1],
							support[2] + bind[2]+1)) + bind[1] + support[1] - 1
			if(!missing(support) && length(support) > 1){
				supLen <- which.max(IRanges::window(bindSpline$y, 2*bindPeak+support[1],
								2*bindPeak+support[2]+1)) +	support[1] - 1 
				bindLen <- bindPeak - 2*supLen
				supPeak <- 2*bindPeak + supLen
				if(verbose) message("Estimated length of support region: ", supLen)
			}
			
			if(verbose) message("Estimated length of binding site: ", bindLen)
		}
	} 
		
	result <- c(bind=bindLen, support=supLen)
		
	if(plot){
		plot(0:max.lag, bindCor, xlab="Lag", ylab="Cross-correlation")
		lines(bindSpline, col=2)
		if(length(bind) > 1) abline(v=bind + support, lty=2)
		if(!missing(support) && length(support) > 1) 
			abline(v=if(length(bind)>1) 2*bindPeak+support else bindLen+2*support, lty=2)
		if(length(bind) > 1) abline(v=bindPeak, lty=3, col=2)
		if(!missing(support) && length(support) > 1) abline(v=supPeak, lty=3, col=3)
	}
	
	result
}


## call nucleosomes. This is just for convenience  
simpleNucCall <- function(data, bind=128, support=17, background=2000, chrLen, ...){
	if(!missing(chrLen))
		callBindingSites(data, bind=bind, support=support, background=background, chrLen=chrLen, ...)
	else
		callBindingSites(data, bind=bind, support=support, background=background, ...)
}

## x is position, y is score. The higher scoring of two overlapping entries is retained
## if y is missing the first entry is chosen
.noOverlap <- function(x, y, minDist){
	result <- numeric(length(x))
	result[1] <- x[1]
	count <- 1
	idx <- 1
	if(length(x) > 1) 
		for(i in 2:length(x)){
			if(x[i] >= result[count] + minDist){
				result[count+1] <- x[i]
				count <- count + 1
				idx <- i
			}
			else if(!missing(y)){
				if(y[i] > y[idx]) result[count] <- x[i]
			}
		}
	length(result) <- count
	result
}

## get sequences for predicted binding sites, removing overlaps by default
exportBindSequence <- function(prediction, reference, bind, overlap=FALSE, file=""){
	## check arguments
	stopifnot(is(prediction, "BindScore"))
	stopifnot(is(reference, "XStringSet"))
	if(missing(bind)) bind <- prediction@functionCall$bind
	
	## predicted binding sites
	peaks <- as.data.frame(prediction)
	peaks <- lapply(levels(peaks$chromosome), function(chr) subset(peaks, chromosome==chr)[,2:3])
	## remove overlap
	if(!overlap)
		peaks <- lapply(peaks, function(x) .noOverlap(x[,1], x[, 2], minDist=bind))
	else peaks <- lapply(peaks, "[[", 1)
	
	start <- lapply(peaks, "-", floor(bind/2))
	end <- lapply(peaks, "+", floor(bind/2))
	
	seqs <- mapply(function(chr, s, e) Views(reference[[chr]], start=s, end=e, 
						names=paste("site", chr, 1:length(s), sep="_")), 1:length(names(prediction)), start, end)	
	
	if(file=="") return(seqs)
	
	files <- sapply(names(prediction), function(chr) paste(file, "_", chr, ".fasta", sep=""))
	mapply(function(x, file, format) writeXStringSet(as(x, "XStringSet"), file), seqs, files)
	
	invisible(seqs)
}

allWords <- function(alphabet, n){
	x <- vector(n, mode="list")
	for(i in 1:n) x[[i]] <- alphabet
	df <- expand.grid(x)
	result <- character(nrow(df))
	for(i in 1:n) result <- paste(result, df[[i]], sep="")
	result
}

allGroups <- function(alphabet, n){
	words <- allWords(alphabet, n)
	groups <- lapply(words, I)
	names(groups) <- words
	
	groups
}


.nucFreq <- function(seq, classes=list(AT=c("AA","AT", "TA", "TT"), GC=c("CC", "CG", "GC", "GG"))){
	patternLen <- nchar(classes[[1]][1])
	counts <- matrix(0, nrow=nchar(seq) - patternLen + 1, ncol=length(classes))
	colnames(counts) <- names(classes)
	for(i in 1:(nchar(seq) - patternLen + 1)){
		counts[i,] <- counts[i,] + sapply(classes, function(x) substr(seq, i, i + patternLen - 1) %in% x)
	}
	counts
}

views2nucFreq <- function(x, classes=list(AT=c("AA","AT", "TA", "TT"), GC=c("CC", "CG", "GC", "GG"))){
	stopifnot(is(x, "XStringViews"))
	patternLen <- nchar(classes[[1]][1])
	## TODO: check pattern length
	
	counts <- matrix(0, nrow=max(width(x)) - patternLen + 1, ncol=length(classes))
	colnames(counts) <- names(classes)
	
	for(i in 1:length(x)){
		f <- .nucFreq(toString(x[i]), classes)
		if(nrow(f) == nrow(counts))	counts <- counts + f
		## TODO: not sure how to best handle this. Only use for uniform windows!
	}
	
	counts
}

.twoWayTest <- function(x, classes=c("A", "C", "G", "T"), df, blocked=FALSE){
	if(is.list(classes)){
		classes1 <- classes[[1]]
		classes2 <- classes[[2]]
	}
	else{
		classes1 <- classes2 <- classes
	}
	dfOffset <- 0
	
	x <- matrix(x, nrow=length(classes1))
	rownames(x) <- classes1
	colnames(x) <- classes2
	
	if(blocked > 0){
		nblock <- nrow(x)/4
		#overlap <- 2*nchar(classes1[1]) - nchar(colnames(x)[1])
		
		blockIdx <- numeric(nblock)
		block <- 1
		for(i in seq(1,length(classes1), by=4)){
			pattern <- substr(classes1[i], nchar(classes1[i]) - blocked + 1, nchar(classes1[i]))
			for(j in 1:length(classes2))
				if(substr(classes2[j], 1, blocked) == pattern){
					blockIdx[block] <- j
					break
				}
			block <- block + 1
		}
		
		result <- vector(nblock, mode="list")
		for(i in 1:nblock){
			result[[i]] <- Recall(x[((i-1)*4+1):(i*4), (blockIdx[i]):(blockIdx[i]+3)])
		}
		
		stat <- sum(sapply(result, "[[", "statistic"))
		df <- sum(sapply(result, "[[", "df"))
		comp <- matrix(0, nrow=nrow(x), ncol=ncol(x))
		for(i in 1:nblock)
			comp[((i-1)*4+1):(i*4), (blockIdx[i]):(blockIdx[i]+3)] <- result[[i]]$components
	}
	else{
		relFreq <- x/sum(x)
		
		m <- list()
		m[[1]] <- rowSums(relFreq)
		m[[2]] <- colSums(relFreq)
		expected <- matrix(apply(expand.grid(m),1,prod)*sum(x),nrow=nrow(x))
		isNull <- sapply(expected, function(x) isTRUE(all.equal(x,0)))
		
		comp <- ifelse(isNull, 0, sign(x - expected) * sqrt((x - expected)^2/expected))
		stat <- sum(comp^2)
		dfOffset <- sum(isNull)
		if(dfOffset) warning("Table contains zero counts. Result may be unreliable.")
	}
	if(missing(df) && blocked) df <- 9 * nblock -dfOffset
	if(missing(df)) df <- (nrow(x)-1)*(ncol(x)-1) - dfOffset
	
	pvalue <- pchisq(stat, df, lower.tail=FALSE)
	
	list(statistic=stat, pvalue=pvalue, components=comp, df=df)
}


.oneWayTest <- function(freq, ref){
	if(any(ref > 1)) ref <- ref/sum(ref)
	expected <- ref * sum(freq[1,])
	comp <- t(apply(freq, 1, function(x) sign(x - expected) * sqrt((x - expected)^2/expected)))
	stat <- rowSums(comp^2)
	pvalue <- pchisq(stat, ncol(freq) -1, lower.tail=FALSE )
	
	list(statistic=stat, pvalue=pvalue, components=comp)
}

diNucTest <- function(dinucFreq, reference){
	order <- c("AA","CA","GA", "TA","AC","CC","GC","TC","AG","CG","GG","TG", "AT","CT","GT","TT")
	dinucFreq <- dinucFreq[, order]
	
	if(missing(reference)){
		## two-way test
		#test <- apply(dinucFreq, 1, .diNucChisqTest2)
		test <- apply(dinucFreq, 1, .twoWayTest)
		stat <- sapply(test, "[[", "statistic")
		pvalue <- sapply(test, "[[", "pvalue")
		comp <- t(sapply(test, function(x) c(x[["components"]])))
		colnames(comp) <- order
	}
	else{
		## one-way test
		reference <- reference[order]
		test <- .oneWayTest(dinucFreq, reference)
		comp <- test$components
		stat <- test$statistic
		pvalue <- test$pvalue
		colnames(comp) <- order
	}
	
	list(statistic=stat, pvalue=pvalue, components=comp)
}

.triNucStat <- function(nuc, mononucFreq, dinucFreq, trinucFreq, n){
	expected <- dinucFreq[-1, substr(nuc,2,3)]/mononucFreq[-c(1,nrow(mononucFreq)), substr(nuc,2,2)]
	observed <- trinucFreq[ , nuc]/n
	
	(observed - expected)^2/expected
}
triNucTest2 <- function(trinucFreq){
	order <- mkAllStrings(c("A", "C", "G", "T"), 3)
	trinucFreq <- trinucFreq[, order]
	subStr1 <- mkAllStrings(c("A", "C", "G", "T"), 2, "left")
	subStr2 <- mkAllStrings(c("A", "C", "G", "T"), 2, "right")
	classes <- expand.grid(subStr1, subStr2)
	
	table <- matrix(0, ncol=nrow(classes), nrow=nrow(trinucFreq))
	idx <- numeric(ncol(trinucFreq))
	names(idx) <- order
	for(nuc in order){
		idx[nuc] <- which(classes[,1] == substr(nuc, 1, 2) & classes[,2] == substr(nuc, 2, 3))
		table[,idx[nuc]] <- trinucFreq[, nuc]
	}
	## df = 36 = 4*9 (9 df from each of the 4 blocks)
	test <- apply(table, 1, .twoWayTest, classes=list(subStr1,subStr2), df=36, blocked=TRUE) 
	
	## restructure result
	stat <- sapply(test, "[[", "statistic")
	pval <- sapply(test, "[[", "pvalue")
	comp <- t(sapply(test, function(x) x$components[idx]))
	colnames(comp) <- order
	
	list(statistic=stat, pvalue=pval, components=comp)
}

triNucTest1 <- function(trinucFreq, reference){
	## one-way test
	reference <- reference[colnames(trinucFreq)]
	.oneWayTest(trinucFreq, reference)
}

## two way test for arbitrarily long nucleotides
oligoNucTest <- function(oligoFreq, alphabet=c("A", "C", "G", "T"), subLen){
	if(missing(subLen)) subLen <- nchar(colnames(oligoFreq)[1]) - 1
	oligo <- nchar(colnames(oligoFreq)[1])
	
	## row and column classes
	subStr1 <- mkAllStrings(alphabet, subLen, "left")
	subStr2 <- mkAllStrings(alphabet, subLen, "right")
	classes <- expand.grid(subStr1, subStr2)
	
	## expand table, each row contains all entries of table for one position
	table <- matrix(0, ncol=nrow(classes), nrow=nrow(oligoFreq))
	idx <- numeric(ncol(oligoFreq))
	names(idx) <- colnames(oligoFreq)
	for(nuc in colnames(oligoFreq)){
		idx[nuc] <- which(classes[,1] == substr(nuc, 1, subLen) & classes[,2] == substr(nuc, oligo-subLen+1, oligo))
		table[,idx[nuc]] <- oligoFreq[, nuc]
	}
	
	## position specific two-way chi^2 test
	test <- apply(table, 1, .twoWayTest, classes=list(subStr1,subStr2), blocked=2*subLen - oligo)

	## restructure result
	stat <- sapply(test, "[[", "statistic")
	pval <- sapply(test, "[[", "pvalue")
	comp <- t(sapply(test, function(x) x$components[idx]))
	colnames(comp) <- colnames(oligoFreq)
	df <- sapply(test, "[[", "df")
	
	list(statistic=stat, pvalue=pval, components=comp, df=df)
}

plotFreq <- function(data, nameLen=1, alphabet=c("A", "C", "G", "T"), 
		xlim=c(-200, 200), bind=c(-73, 73), support=10, range=c(-500,499)){
	par.old <- par(no.readonly=TRUE)
	layout.mat <- matrix(1:16,ncol=4)
	par(mar=rep(0, 4))
	par(oma=rep(4.1, 4))
	layout(layout.mat)
	
	rowNames <- mkAllStrings(alphabet, nameLen, "left")
	colNames <- mkAllStrings(alphabet, nameLen, "right")
	
#	rowNames <- unique(sapply(colnames(data), function(x) substr(x, 1, nameLen)))
#	colNames <- unique(sapply(colnames(data), function(x) substr(x, nchar(x) - nameLen + 1, nchar(x))))

	bindLen <- bind[2]-bind[1]
	ticks <- c(seq(sign(bind[1])*ceiling(bindLen/5)*5, xlim[1], by=-ceiling(bindLen/5)*5), 
			bind[1], (bindLen)/2 + bind[1], bind[2], 
			seq(sign(bind[2])*ceiling(bindLen/5)*5, xlim[2], by=ceiling(bindLen/5)*5))
	
	## determine row and column indices for plotting of 4x4 blocks
	row <- numeric()
	for(i in seq(4, length(rowNames), by=4)) row <- c(row, rep((i-3):i, times=4))
	row <- rep(row, times=ceiling(length(colNames)/4))
	
	col <- numeric()
	for(i in seq(4, length(colNames), by=4)) col <- c(col, rep(rep((i-3):i, each=4), 
						times=ceiling(length(colNames)/4)))
	
	## if plotting a blocked design, i.e. row and column names overlap, remove all
	## non matching combinations (structural 0s)
	overlap <- 2*nameLen - nchar(colnames(data)[1])
	if(overlap){
		match <- mapply(function(r, c) 
					substr(rowNames[r], nameLen-overlap +1, nameLen) == substr(colNames[c], 1, overlap), 
				row, col)
		row <- row[match]
		col <- col[match]
	}
	
	idx <- 1:length(row)
	for(i in idx){
		nuc <- paste(rowNames[row[i]], substr(colNames[col[i]], overlap + 1, nameLen), sep="")
		
		plot(range[1]:range[2], data[,nuc], ylim=range(data),xlim=xlim,xaxt="n", yaxt="n", xlab="", ylab="")
		lines(bind, c(0,0), col=rgb(0,1,0, 0.75), lwd=3, lend=1)
		lines(bind + c(-support, support), c(0,0), col=rgb(0,1,0, 0.75), lwd=1, lend=1)
		sp <- smooth.spline(range[1]:range[2], data[,nuc], nknots=200)
		lines(sp, col=2)
		
		if(col[i]%%4 == 1){
#			if(length(rowNames)*length(colNames) == ncol(data))
			mtext(rowNames[row[i]], 2, 2)
#			else ## blocked layout
#				mtext(rowNames[floor(i/16)*4 + ((i-1)%%16)+1], 2, 2)
		}
		if(col[i] %% 4 == 0){
			axis(4)
		}
		if(row[i] %% 4 == 1){
#				if(length(rowNames)*length(colNames) == ncol(data))
			mtext(colNames[col[i]], 3, 2)
#				mtext(colNames[floor(i/16) * 4 + (((((i - 1)%%16) + 1)/4) + 1)], 3, 2)
		}
		if(row[i] %% 4 == 0){
			axis(1, at=ticks)
		}
		if(row[i]*col[i] %% 16 == 0){
			mtext("Distance from nucleosome centre", 1, 3, outer=TRUE)
			mtext("Score", 4, 3, outer=TRUE)
		}
	}
	
#	for(i in 1:ncol(data)){
#		plot(range[1]:range[2], data[,i], ylim=range(data),xlim=xlim,xaxt="n", yaxt="n", xlab="", ylab="")
#		lines(bind, c(0,0), col=rgb(0,1,0, 0.75), lwd=3, lend=1)
#		lines(bind + c(-support, support), c(0,0), col=rgb(0,1,0, 0.75), lwd=1, lend=1)
#		sp <- smooth.spline(range[1]:range[2], data[,i], nknots=200)
#		lines(sp, col=2)
#		row <- (i-1)%%length(rowNames) + 1
#		col <- ceiling(i/length(rowNames))
#		if(col%%4 == 1){
##			if(length(rowNames)*length(colNames) == ncol(data))
#				mtext(rowNames[row], 2, 2)
##			else ## blocked layout
##				mtext(rowNames[floor(i/16)*4 + ((i-1)%%16)+1], 2, 2)
#		}
#		if((((i-1)%%16)+1) %in% 13:16){
#			axis(4)
#		}
#		if((((i-1)%%16)+1) %in% c(1,5,9,13)){
#			if(length(rowNames)*length(colNames) == ncol(data))
#				mtext(colNames[floor(i/16)+1], 3, 2)
#			mtext(colNames[floor(i/16) * 4 + (((((i - 1)%%16) + 1)/4) + 1)], 3, 2)
#		}
#		if((((i-1)%%16)+1) %in% c(4,8,12,16)){
#			axis(1, at=ticks)
#		}
#		if(i %% 16 == 0){
#			mtext("Distance from nucleosome centre", 1, 3, outer=TRUE)
#			mtext("Score", 4, 3, outer=TRUE)
#		}
#	}
	
	par(par.old)
}

.oligonucSim <- function(freq1, freq2, oligo=3, iter=3000){
	freq1 <- freq1[mkAllStrings(c("A", "C", "G", "T"), oligo-1, "left")]
	freq2 <- freq2[mkAllStrings(c("A", "C", "G", "T"), oligo-1, "right")]
	
	## calculate oligonucleotide frequencies
	oligonuc <- mkAllStrings(c("A", "C", "G", "T"), oligo)
	oligonucFreq <- numeric(length(oligonuc))
	names(oligonucFreq) <- oligonuc
	for(nuc in oligonuc){
		oligonucFreq[nuc] <- freq1[substr(nuc, 1, oligo-1)] * freq2[substr(nuc, 2, oligo)]
	}
	
	## generate oligonucleotide counts
	oligonucSample <- sample(oligonuc, iter, prob=oligonucFreq, replace=TRUE)
	oligonucSample <- table(oligonucSample)
	
	## arrange oligonucleotides into blocks
	classes <- expand.grid(names(freq1), names(freq2))
	table <- numeric(nrow(classes))
	idx <- numeric(length(oligonuc))
	names(idx) <- oligonuc
	for(nuc in oligonuc){
		idx[nuc] <- which(classes[,1] == substr(nuc, 1, oligo-1) & classes[,2] == substr(nuc, 2, oligo))
		table[idx[nuc]] <- ifelse(is.na(oligonucSample[nuc]), 0, oligonucSample[nuc])
	}
		
	test <- .twoWayTest(table, classes=list(names(freq1), names(freq2)), blocked=oligo-2)
	test$statistic
}

## print head and tail of a vector
.headTail <- function(x, n = 6L){
	start <- paste(head(x, n), collapse = ", ")
	middle <- if(length(x) == 2*n + 1) {
		paste(", ", x[n+1], ", ", sep="") 
	}else if(length(x) > 2*n + 1) ", ..., " else if(length(x) < 2*n + 1 && length(x) > n) ", "
	end <- if(length(x) > n) paste(tail(x, min(n, length(x) - n)), collapse = ", ")
	paste(start, middle, end, sep = "")
} 

## helper function to print common information about objects of classes defined in this package
.showHeader <- function(object){
	cat("Class:", class(object),"\n")
	if("functionCall" %in% slotNames(object)){
		cat("Generated by: ")
		show(object@functionCall)
	}
	cat("Chromosomes: ", length(chrLength(object)), " (", .headTail(names(object), 2) ,")\n", sep="")
	cat("Genome length:", sum(chrLength(object)), "\n")
}

## calculate p-values
.score2pval <- function(score, null){
	
	scores <- c(lapply(score, function(y) 
						as.numeric(y[!is.na(y)])), recursive=TRUE)
	pval <- pnorm(scores, mean=null["mean"], sd=null["sd"], lower.tail=FALSE)
	pval <- p.adjust(pval, "BH")
	## split back into chromosomes
	len <- sapply(score, function(y) sum(!is.na(y)))
	end <- cumsum(len)
	start <- c(1, end[-length(end)]+1)
	
	result <- vector(mode="list", length(score))
	for(i in 1:length(result)){
		result[[i]][!is.na(score[[i]])] <- pval[start[i]:end[i]]
	}
	result
}

## append trailing zeros to get full chromosome length
## returns RleList with read counts for both strands
.fixCounts <- function(counts, totalLen){
	countLen <- sapply(counts, length)
	if(countLen[1] < totalLen){
		counts[[1]]@lengths <- c(runLength(counts[[1]]), as.integer(totalLen) - countLen[1])
		counts[[1]]@values <- c(runValue(counts[[1]]), 0L)
	}
	if(countLen[2] < totalLen){
		counts[[2]]@lengths <- c(runLength(counts[[2]]), as.integer(totalLen) - countLen[2])
		counts[[2]]@values <- c(runValue(counts[[2]]), 0L)
	}
	
	RleList(counts[[1]], counts[[2]])
}
