% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_match_class.R
\name{summary,graphMatch-method}
\alias{summary,graphMatch-method}
\title{Summary methods for graphMatch objects}
\usage{
\S4method{summary}{graphMatch}(object, A = NULL, B = NULL, true_label = NULL, directed = NULL)
}
\arguments{
\item{object}{graphMatch object}

\item{A}{igraph or matrix-like object}

\item{B}{igraph or matrix-like object}

\item{true_label}{the true correspondence (if available)}

\item{directed}{whether to treat the graphs as directed (TRUE) or not
directed (FALSE) default is NULL which will treat the graphs as directed if
either adjacency matrix is not symmetric.}
}
\value{
\code{summary} returns the graph matching formula, and a summary of
  graph matching results including the number of matches, the number of
  correct matches (if the true correspondence is available), and common
  edges, missing edges, extra edges, common non-edges and the objective
  function value.
}
\description{
Summary methods for graphMatch objects
}
\examples{
set.seed(123)
graphs <- sample_correlated_gnp_pair(20, .9, .3)
A <- graphs$graph1
B <- graphs$graph2
match <- gm(A, B, 1:4, method = "percolation")

summary(match, A, B)
summary(match, A, B, true_label = 1:20) # also output the number of correct matches

}
