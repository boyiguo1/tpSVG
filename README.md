
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tpSVG <img src="logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R build
status](https://github.com/boyiguo1/tpSVG/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/boyiguo1/tpSVG/actions)
<!-- badges: end -->

The goal of `tpSVG` is to detect and visualize spatial variation in the
gene expression for spatially resolved transcriptomics data analysis.
Specifically, `tpSVG` introduces a family of count-based models, with
generalizable parametric assumptions such as Poisson distribution or
negative binomial distribution. In addition, comparing to
crmarkdown::pandoc_version()urrently available count-based model for
spatially resolved data analysis, the `tpSVG` models improves
computational time, and hence greatly improves the applicability of
count-based models in SRT data analysis.

## Installation

### GitHub

You can install the development version of tpSVG from
[GitHub](https://github.com/boyiguo1/tpSVG) with:

``` r
#' Install devtools package if not already installed
if (required(devtools)) install.packages(package_name)
devtools::install_github("boyiguo1/tpSVG")
```

If you have R version before v4.4 and would like to install tpSVG, you
can follow

    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("boyiguo1/tpSVG@pre-R4.4")

> WARNING: The purpose of having the branch pre-R4.4 is to allow users
> to use escheR before the formal release of R 4.4 and during the early
> stage of R 4.4 release. This branch will not be update with any
> further development beyond escheR v0.99.1. We recommend users to
> update their R versions up to date.

### Bioconductor (pending)

The package is currently submitted to Bioconductor for
[review](https://github.com/Bioconductor/Contributions/issues/3264).
Once the package is accepted by Bioconductor, you can install the latest
release version of `tpSVG` from Bioconductor via the following code.
Additional details are shown on the Bioconductor page.

``` r
# NOTE: The package is under-review with bioconductor.
#       The following code section will work once the package is accepted.

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("tpSVG")
```

The latest development version can also be installed from the `devel`
version of Bioconductor or from GitHub following

``` r
BiocManager::install(version = "devel")
```

## Tutorial

Please find an end-to-end tutorial at
<https://boyi-guo.com/tpSVG/articles/intro_to_tpSVG.html>.

## Frequently asked questions

**Implementation Questions**

- What are the data structures that `tpSVG` current supports?

  *As of `tpSVG v0.99.1`, the data structure `tpSVG` supports includes
  [`SpatialExperiments`](https://bioconductor.org/packages/release/bioc/html/SpatialFeatureExperiment.html)
  (and packages extending `SpatialExperiments`,
  e.g. [`SpatialFeatureExperiments`](https://bioconductor.org/packages/release/bioc/html/SpatialFeatureExperiment.html))
  and `data.frame`. Please find example via
  [supported_data_structure](https://boyi-guo.com/escheR/articles/supported_data_structure.html).
  Due to limited resources, we regret that we won’t provides direct
  accessibility to other pipelines, e.g. `suerat`.*

- What types of spatially-resolved transcriptomics (SRT) data that
  `tpSVG` supports?

  *Both sequenced-based SRT and image-based SRT data are supported by
  `tpSVG`. For more details, please refer to the vignette
  \[supported_data_structure\]\](<https://boyi-guo.com/tpSVG/articles/supported_data_structure.html#image-based-srt-in-spatialexperiment-e-g--spatialfeatureexperiment>).*

- Can I use other scale factor as offset in the count-model?

  *Yes, just remember to take log for the offset term. In the vignettes,
  the offset of the model is default to library size, i.e. the total
  number of molecular in a spot/cell, but the count models should be
  compatible to other definition of scale factor in theory.*

**Theoretical Questions**

- What is the difference between modeling log transformed data and count
  data?

  *Count data is the natural form of gene expression data when it is
  collected and quantified. While log-transformation providess shortcuts
  to model (normalized) count data using well-studied Gaussian
  distribution, it distorts the lowly expressed gene and causes analytic
  biases.*
