![](http://bioconductor.org/shields/availability/devel/ChIPseqR.svg)
![](http://bioconductor.org/shields/build/devel/bioc/ChIPseqR.svg)
![](http://bioconductor.org/shields/years-in-bioc/ChIPseqR.svg)
![](http://bioconductor.org/shields/downloads/ChIPseqR.svg)

# ChIPseqR
ChIPseqR identifies protein binding sites from ChIP-seq and nucleosome positioning experiments. 
The model used to describe binding events was developed to locate nucleosomes but should flexible enough to 
handle other types of experiments as well.

# Installation
Both the [development](http://bioconductor.org/packages/devel/bioc/html/ChIPseqR.html) and 
[release](http://bioconductor.org/packages/release/bioc/html/ChIPseqR.html) version of this R package
are available through [Bioconductor](http://bioconductor.org/). Use the *BiocManager* to install the
package and its dependencies from within R.

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ChIPseqR")
```
