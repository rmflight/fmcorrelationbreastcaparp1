[![Build Status](https://travis-ci.org/rmflight/fmcorrelationbreastcaparp1.svg?branch=master)](https://travis-ci.org/rmflight/fmcorrelationbreastcaparp1)

# Fondufe-Mittendorf Breast Cancer PARP1 Correlation Analysis Package

This package contains functions and vignettes used in generating correlation results for the publication:

**Genome-wide profiling of PARP1 reveals an interplay with gene regulatory regions and DNA methylation**, Nalabothula Narasimharao, Taha Al-jumaily, Robert M. Flight, Shao Xiaorong, Hunter N. B. Moseley, F. Lisa Barcellos, Yvonne Fondufe-Mittendorf, Submitted

## Contents

* Functions used in the analysis
* Vignettes providing the full calculations of the correlation of nucleosome associated PARP1 reads in the MCF-7 and MDA-MB123 cell lines with
  * [CTCF ChIP-Seq](https://rmflight.github.io/fmcorrelationbreastcaparp1/parp1_ctcf.html)
  * [Methylation reads](https://rmflight.github.io/fmcorrelationbreastcaparp1/parp1_methylation.html)
  * [Histone Marks ChIP-Seq](https://rmflight.github.io/fmcorrelationbreastcaparp1/parp1_histone_marks.html)
  * [Transcript Expression](https://rmflight.github.io/fmcorrelationbreastcaparp1/parp1_expression.html)
  
In addition, the correlations are stored in plain text files available [here](https://github.com/rmflight/fmcorrelationbreastcaparp1/tree/master/inst/correlation_tables). These values were reported in the paper.

## Citation

If this package is used for other analyses, it should be cited as:

**fmcorrelationbreastcaparp1: Functions for calculating the correlations in PARP1 nucleosome data in breast cancer cell lines**, R. M. Flight, Y. Fondufe-Mittendorf, H. N. B. Moseley doi:10.6084/m9.figshare.1267544

## Installation

To install this package, you should use `devtools`:

```
library("devtools")
install_github("rmflight/fmcorrelationbreastcaparp1")
```

## Recreating Vignettes

If you want to completely repeat the analysis in the vignettes, you will also need to install the [data package](https://github.com/rmflight/fmdatabreastcaparp1), and you will want to clone this package locally and run `devtools::build_vignettes`.

On the command line, first install the data package, and clone this analysis package

```
# clone the packages
git clone https://github.com/rmflight/fmdatabreastcaparp1.git
git clone https://github.com/rmflight/fmcorrelationbreastcaparp1.git

# create the data directory
mkdir fmdatabreastcaparp1/data

# pull down the data
wget http://downloads.figshare.com/article/public/1266451.zip
unzip 1266451.zip -d fmdatabreastcaparp1/data

# start R and use devtools to install the packages
R
library(devtools)
install("fmdatabreastcaparp1")
install("fmcorrelationbreastcaparp1")
```

Then in `R`, use `devtools` to re-build the vignettes:

```
# assumes you started R in the cloned directory
devtools::build_vignettes()
```

This will generate the correlation files as well, hopefully they will give the same results as are already stored.

## Memory

Re-generating the vignettes is a memory hog. Tracking memory usage during `build_vignettes` resulted in ~16GB used at peak areas. Keep that in mind if you want to redo the vignettes.
