#' @title Find the largest common connected subgraph (LCCS) of two graphs
#'
#' @description Find the largest common connected subgraphs of
#' two matched graphs, which is an induced connected subgraph of both graphs 
#' that has as many vertices as possible.
#' The \code{largest_cc} function returns the largest connected subgraph of a single graph.
#'
#' @param A A matrix or an 'igraph' object. See \link{check_graph}. Must be single-layer.
#' @param B A matrix or an 'igraph' object. See \link{check_graph}. Must be single-layer.
#' @param min_degree A number. Defines the level of connectedness of the 
#' obtained largest common connected subgraph. The induced subgraph is 
#' a graph with a minimum vertex-degree of at least min_degree.
#'
#' @rdname largest_common_cc
#'
#' @return \code{largest_common_cc} returns the common largest connected subgraphs of
#' two aligned graphs in the 'igraph' object form and a logical vector indicating which vertices in
#' the original graphs remain in the induced subgraph.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.7, p =  0.2)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # put no constraint on the minimum degree of the common largest conncect subgraph
#' lccs1 <- largest_common_cc(g1, g2, min_degree = 1)
#' # induced subgraph
#' lccs1$g1
#' lccs1$g2
#' # label of vertices of the induced subgraph in the original graph
#' igraph::V(g1)[lccs1$keep]
#'
#' # obtain a common largest connect subgraph with each vertex having a minimum degree of 3
#' lccs3 <- largest_common_cc(g1, g2, min_degree = 3)
#' @export
#'
largest_common_cc <- function(A, B, min_degree = 1){
  graph_pair <- check_graph(A, B, as_list = FALSE)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  A <- igraph::graph_from_adjacency_matrix(as.matrix(A))
  B <- igraph::graph_from_adjacency_matrix(as.matrix(B))
  keep <- rep(TRUE, igraph::vcount(A))
  while (!(igraph::is_connected(A) && igraph::is_connected(B))){
    cc1 <- igraph::components(A)
    cc2 <- igraph::components(B)

    lcc1 <- which.max(cc1$csize)
    lcc2 <- which.max(cc2$csize)

    vlcc1 <- cc1$membership == lcc1
    vlcc2 <- cc2$membership == lcc2

    vlcc <- vlcc1 & vlcc2
    keep[keep] <- vlcc

    A <- igraph::induced_subgraph(A, igraph::V(A)[vlcc])
    B <- igraph::induced_subgraph(B, igraph::V(B)[vlcc])

    if(min_degree > 1){
      good_deg <- (igraph::degree(A) >= min_degree) &
        (igraph::degree(B) >= min_degree)
      keep[keep] <- good_deg

      A <- igraph::induced_subgraph(A, good_deg)
      B <- igraph::induced_subgraph(B, good_deg)
    }
  }

  list(g1 = A, g2 = B, keep = keep)
}

#' @rdname largest_common_cc
#' @export
largest_cc <- function(A){
  g <- igraph::graph_from_adjacency_matrix(as.matrix(check_single_graph(A)))
  c <- igraph::components(g)
  lc <- which.max(c$csize)
  keep <- c$membership == lc
  list(g = induced_subgraph(g, V(g)[keep]), keep = keep)
}
