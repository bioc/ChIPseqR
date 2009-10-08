# New generic functions
# 
# Author: Peter Humburg
###############################################################################

## create strand specific read counts
setGeneric("strandPileup", def=function(aligned, chrLen, ...) standardGeneric("strandPileup"))

## call binding sites
## Convenience function for one step analysis
## This function combines the different steps required to identify nucleosomes into a single function.
## 'data' is either an object with mapped reads or a list with read counts (as returned by strandPileup).
## Defaults provided for the various parameters should be reasonable for nucleosome data.
## An attempt is made to determine the length of binding and support regions from the data if more than one value
## is provided, see getBindLen for details.
## Returns a list with components 
##    binding: a data.frame with predicted binding site centres and associated scores and p-values
##    score: a list with all calculated scores (one numeric vector per chromosome)
##    pval: a list with all calculated p-values (one numeric vector per chromosome)
setGeneric("callBindingSites", def=function(data, ...) standardGeneric("callBindingSites"))


## accessors for slots of BindScore
setGeneric("pvalue", def=function(x, ...) standardGeneric("pvalue"), useAsDefault=function(x, ...) x@pvalue)
#setGeneric("score", def=function(x, ...) standardGeneric("score"), useAsDefault=function(x, ...) x@score)
setGeneric("peaks", def=function(x, ...) standardGeneric("peaks"), useAsDefault=function(x, ...) x@peaks)
setGeneric("binding", def=function(x, ...) standardGeneric("binding"), useAsDefault=function(x, ...) x@binding)
setGeneric("support", def=function(x, ...) standardGeneric("support"), useAsDefault=function(x, ...) x@support)
setGeneric("cutoff", def=function(x, ...) standardGeneric("cutoff"), useAsDefault=function(x, ...) x@cutoff)
setGeneric("cutoff<-", def=function(x, ..., value) standardGeneric("cutoff<-"), 
		useAsDefault=function(x, ...,  value) x@cutoff <- value)
setGeneric("nullDist", def=function(x, ...) standardGeneric("nullDist"), useAsDefault=function(x, ...) x@nullDist)
setGeneric("nullDist<-", def=function(x, ..., value) standardGeneric("nullDist<-"), 
		useAsDefault=function(x, ..., value) x@nullDist <- value)

## explicit conversion between compressed and expanded representations
setGeneric("compress", def=function(x, ...) standardGeneric("compress"))
setGeneric("decompress", def=function(x, ...) standardGeneric("decompress"), useAsDefault=function(x, ...) x)

## create S4 generics for some existing (S3) methods
setGeneric("lapply")
setGeneric("sapply")
setGeneric("as.data.frame")
setGeneric("head", useAsDefault=utils::head)
setGeneric("tail", useAsDefault=utils::tail)
setGeneric("plot", useAsDefault=graphics::plot)

## summary functions
##    number of reads in sample
setGeneric("nreads", def=function(x, ...) standardGeneric("nreads"))
##    length of chromosomes
setGeneric("chrLength", def=function(x, ...) standardGeneric("chrLength"))