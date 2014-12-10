# Fondufe-Mittendorf Breast Cancer PARP1 Correlation Analysis Package

This package contains functions and vignettes used in generating correlation results for the publication:

**Genome-wide profiling of PARP1 reveals an interplay with gene regulatory regions and DNA methylation**, Nalabothula Narasimharao, Taha Al-jumaily, Robert M. Flight, Shao Xiaorong, Hunter N. B. Moseley, F. Lisa Barcellos, Yvonne Fondufe-Mittendorf, Submitted

## Contents

* Functions used in the analysis
* Vignettes providing the full calculations of the correlation of nucleosome associated PARP1 reads in the MCF-7 and MDA-MB123 cell lines with
  * [CTCF ChIP-Seq](ctcf_vignette)
  * [Methylation reads](methylation_vignette)
  * [Histone Marks ChIP-Seq](histone_vignette)
  * [Transcript Expression](expression_vignette)
  
In addition, the correlations are stored in plain text files available [here](https://github.com/rmflight/fmanalysisbreastcaparp1/tree/master/inst/correlation_tables). These values were reported in the paper.

## Citation

If this package is used for other analyses, it should be cited as:

**fmanalysisbreastcaparp1: Functions for analyzing PARP1 nucleosome data on breast cancer cell lines**, R. M. Flight, Y. Fondufe-Mittendorf, H. N. B. Moseley doi:xxxxxx

## Installation

To install this package, you should use `devtools`:

```
library("devtools")
install_github("rmflight/fmanalysisbreastcaparp1")
```

## Recreating Vignettes

If you want to completely repeat the analysis in the vignettes, you will also need to install the [data package](datalink), and you will want to clone this package locally and run `devtools::build_vignettes`.

On the command line, first clone this analysis package:

```
git clone https://github.com/rmflight/fmanalysisbreastcaparp1.git
```

Then in `R`, install the data package, and re-build the vignettes:

```
# assumes you started R in the cloned directory
library("devtools")
install_github("rmflight/fmdatabreastcaparp1")
library("fmdatabreastcaparp1")
build_vignettes()
```

This will generate the correlation files as well, hopefully they will give the same results as are already stored.
