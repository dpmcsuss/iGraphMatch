# iGraphMatch

iGraphMatch is a R package for graph matching. The package works for both igraph objects and matrix objects. You provide the adjacency matrices of two graphs and some other information you might know, choose the graph matching method, and it returns the graph matching results. iGraphMatch also provides a bunch of useful functions you might need during the process of graph matching.

Installation
------------

This is the development version of the package. This can be installed via devtools with
``` r
# install.packages("devtools")
devtools::install_github("dpmcsuss/iGraphMatch", ref = "dev")
```

You might want to try the `dev` branch since it might be more stable at the moment.

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
