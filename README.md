<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/dpmcsuss/iGraphMatch/branch/dev/graph/badge.svg)](https://codecov.io/gh/dpmcsuss/iGraphMatch?branch=dev)
[![CRAN status](https://www.r-pkg.org/badges/version/iGraphMatch)](https://CRAN.R-project.org/package=iGraphMatch)
[![R-CMD-check](https://github.com/dpmcsuss/iGraphMatch/workflows/R-CMD-check/badge.svg)](https://github.com/dpmcsuss/iGraphMatch/actions)
<!-- badges: end -->

<!-- [![Build Status](https://travis-ci.com/dpmcsuss/iGraphMatch.svg?branch=dev)](https://travis-ci.com/dpmcsuss/iGraphMatch) -->

# iGraphMatch

iGraphMatch is a R package for graph matching. The package works for both igraph objects and matrix objects. You provide the adjacency matrices of two graphs and some other information you might know, choose the graph matching method, and it returns the graph matching results. iGraphMatch also provides a bunch of useful functions you might need during the process of graph matching.

Installation
------------
``` r
# install.packages("devtools")
devtools::install_github("dpmcsuss/iGraphMatch")
```

This package is still in development. The `dev` branch of the package might be more stable. This can be installed via devtools with

``` r
# install.packages("devtools")
devtools::install_github("dpmcsuss/iGraphMatch", ref = "dev")
```

Usage
------------

After installation, we recommend loading `igraph` first.
``` r
library(igraph)
library(iGraphMatch)
``` 

Documentation
------------

Documentation can be found via the help in R and at https://dpmcsuss.github.io/iGraphMatch/.



Author and Contributors
-----------------------

The primary authors for the package are Vince Lyzinski, Zihuan Qiao, and Daniel Sussman.

Joshua Agterberg, Lujia Wang, and Yixin Kong also provided important contributions. We also want to thank all of our users, especially Youngser Park, for their feedback and patience as we continue to develop the package.


Support
-------

This work was supported in part by grants from DARPA (FA8750-20-2-1001 and FA8750-18-0035) and from MIT Lincoln Labs.
