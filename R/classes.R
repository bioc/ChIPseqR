# S4 class definitions
# 
# Author: Peter Humburg
###############################################################################


## Class representing strand specific read counts
setClass("ReadCounts", representation=representation(counts="list"), prototype=prototype(counts=list()))

## Class representing scored binding sites
## contains all computed scores and associated p-values as well as the location of predicted binding sites
setClass("BindScore", 
		representation=representation(functionCall="call", score="list", 
				pvalue="list", peaks="list", cutoff="numeric", nullDist="numeric", start="integer"),
		prototype=prototype(score=list(), pvalue=list(), peaks=list(), cutoff=numeric(), 
				nullDist=numeric(), start=0L))

## variants of the above using run-length encoding to reduce memory requirements
setClass("RLEReadCounts", contains="ReadCounts")
setClass("RLEBindScore", contains="BindScore")
