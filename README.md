# iGraphMatch

Overview
--------
iGraphMatch is a R package for graph matching. The package works for both igraph objects and matrix objects. You provide the adjacency matrices of two graphs and some other information you might know, choose the graph matching method, and it returns the graph matching results. iGraphMatch also provides a bunch of useful functions you might need during the process of graph matching.

Installation
------------

Currently only the development version of this package is available. This can be installed via devtools with
``` r
# install.packages("devtools")
devtools::install_github("dpmcsuss/iGraphMatch", auth_token = "AUTH_TOKEN")
```
You will need to generate an authorization token using the following instructions.

> To install from a private repo, generate a personal access token (PAT) in https://github.com/settings/tokens and supply to this argument. This is safer than using a password because you can easily delete a PAT without affecting any others. Defaults to the GITHUB_PAT environment variable.
