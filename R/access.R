# Access and replacement functions for S4 classes
# 
# Author: Peter Humburg
###############################################################################

############################### ReadCounts ####################################

## access
setMethod("[[", c(x="ReadCounts", j="missing"),
		function(x, i, exact=TRUE){
			x@counts[[i, exact=exact]]
		}
)

setMethod("[", c(x="ReadCounts", j="missing", drop="missing"),
		function(x, i){
			x@counts[i]
		}
)

setMethod("$", "ReadCounts",
		function(x, name){
			x@counts[[name, exact=FALSE]]
		}
)

setMethod("names", "ReadCounts",
		function(x){
			names(x@counts)
		}
)

## replacement
setMethod("[[<-", c(x="ReadCounts", j="missing"),
		function(x, i, exact=TRUE, value){
			## ensure that value is an integer matrix
			if(!is.matrix(value) || !is.integer(value) || !(ncol(value) == 2)){
				value <- matrix(as.integer(value), ncol=2)
			}
			colnames(value) <- c("+", "-")
			x@counts[[i, exact=exact]] <- value
			
			invisible(x)
		}
)

setMethod("[<-", c(x="ReadCounts", j="missing"),
		function(x, i, exact=TRUE, value){
			for(j in 1:length(i))
				x[[i[[j]], exact=exact]] <- value[[j]]
			
			invisible(x)
		}
)

setMethod("$<-", "ReadCounts",
		function(x, name, value){
			x[[name, exact=FALSE]] <- value
			invisible(x)
		}
)

setMethod("names<-", "ReadCounts", function(x, value){names(x@counts) <- value; invisible(x)})

############################### BindScore ####################################
## subset results to restrict them to individual chromosomes
setMethod("[[", c(x="BindScore", j="missing"),
		function(x, i, exact=TRUE) x[i]
)
## subset results to restrict them to a subset of an individual chromosome
setMethod("[[", c(x="BindScore", j="numeric"),
		function(x, i, j, exact=TRUE){
			if(is.character(i)) names <- i
			else names <- names(x)[i]
			peaks <- peaks(x)[[i]][peaks(x)[[i]] <= max(j) & peaks(x)[[i]] >= min(j)]
			offset <- x@start
			j <- j - offset + 1
			if(max(j) <= 0) j <- numeric()
			else j[which.min(j)] <- max(min(j), 1)
			j <- j[j > 0]
			scores <- if(length(j) > 0) list(as.numeric(window(score(x)[[i]],min(j),max(j)))) 
					else list(new(class(score(x)[[i]])))
			pvalues <- if(length(j) > 0) list(as.numeric(window(pvalue(x)[[i]],min(j),max(j)))) 
					else list(new(class(pvalue(x)[[i]])))
			new(class(x), x@functionCall, scores, pvalues, list(peaks), cutoff=cutoff(x), 
					nullDist=nullDist(x), names=names, start=min(j)+offset - 1)
		}
)

## subset results by chromosome
setMethod("[", c(x="BindScore", j="missing", drop="missing"),
		function(x, i){
			if(is.character(i)) names <- i
			else names <- names(x)[i]
			new(class(x), x@functionCall, score(x)[i], pvalue(x)[i], peaks(x)[i], 
					cutoff=cutoff(x), nullDist=nullDist(x), names=names)
		}
)

setMethod("names", "BindScore", function(x) names(x@score))

setMethod("names<-", "BindScore", function(x, value){
			names(x@score) <- value
			names(x@pvalue) <- value
			names(x@peaks) <- value
			invisible(x)
		}
)

## access slots
setMethod("score", "BindScore", function(x) x@score)
setMethod("cutoff", "BindScore", 
		function(x, type=c("score", "pvalue")){
			select <- pmatch(type, c("score", "pvalue"), nomatch=0)
			if(all(sapply(select, ">", 0)))
				return(x@cutoff[select])
			else return(NULL)
		}
)
setMethod("cutoff<-", "BindScore",
		function(x, type=c("score", "pvalue"), value){
			select <- pmatch(type, c("score", "pvalue"), nomatch=0)
			if(length(select) != 1 || select == 0)
				stop("Argument 'type' has to be either \"score\" or \"pvalue\".")
			
			if(select == 1){
				newScore <- value
				## find matching p-value. We don't have information about null distribution anymore, 
				## infer from data instead
				chrMin <- sapply(score(x), function(y) min(y[y >= value], na.rm=TRUE))
				chr <- which.min(chrMin)
				pv.idx <- which(score(x)[[chr]] == chrMin[chr])[1]
				newPvalue <- pvalue(x)[[chr]][pv.idx]
			}
			else if(select == 2){
				newPvalue <- value
				## find matching score
				chrMax <- sapply(pvalue(x), function(y) max(y[y <= value], na.rm=TRUE))
				chr <- which.max(chrMax)
				score.idx <- which(pvalue(x)[[chr]] == chrMax[chr])[1]
				newScore <- score(x)[[chr]][score.idx]
			}
			
			x@cutoff[1:2] <- c(newScore, newPvalue)
			x@peaks <- lapply(score(x), pickPeak, newScore, offset=support(x) + ceiling(binding(x)/2))
			
			invisible(x)
		}
)
setMethod("nullDist<-", "BindScore",
		function(x, value){
			## check new value
			stopifnot(length(value) == 2)
			x@nullDist[1:2] <- value
			
			## re-evaluate p-values
			x@pvalue <- .score2pval(score(x), nullDist(x))
			
			## update score cut-off and significant peaks
			cutoff(x, type="pvalue") <- cutoff(x, type="pvalue")
			
			invisible(x)
		}
)

## extract parameter values from function call
setMethod("binding", "BindScore", function(x) x@functionCall$bind)
setMethod("support", "BindScore", function(x) x@functionCall$support)


