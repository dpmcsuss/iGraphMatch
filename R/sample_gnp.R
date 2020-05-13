#' @importFrom igraph sample_gnp sample_sbm graph_from_adjacency_matrix %u% %s% is.igraph V

#' @title Sample correlated G(n,p) random graphs
#'
#' @description Sample a pair of correlated G(n,p) random graphs with correlation between
#' two graphs being \code{rho} and edge probability being \code{p}.
#'
#' @param n An integer. Number of total vertices for the sampled graphs.
#' @param corr A number. The target Pearson correlation between the adjacency matrices
#' of the generated graphs. It must be in open (0,1) interval.
#' @param p A number. Edge probability between two vertices. It must be in open
#' (0,1) interval.
#' @param ncore An integer. Number of core vertices.
#' @param permutation A numeric vector,permute second graph. 
#' @param ... Passed to \code{sample_correlated_gnp_pair} and \code{sample_correlated_gnp_pair_w_junk}.
#'
#' @rdname sample_gnp
#' @return \code{sample_correlated_gnp_pair} returns a list of two igraph object, named
#' \code{graph1} and \code{graph2}, which are two graphs whose adjacency matrix entries
#' correlated with \code{rho}.
#' @examples
#' sample_correlated_gnp_pair(50, 0.3, 0.5)
#' @export
#'
sample_correlated_gnp_pair <- function(n, corr, p, permutation=1:n, ...){
  
  # Make the first graph
  graph1 <- sample_gnp(n,p,...)
  
  # Make two graphs which will be used to make the
  # second graph
  Z0 <- sample_gnp(n,p*(1-corr),...)
  Z1 <- sample_gnp(n,p+corr*(1-p),...)
  
  graph2 <- Z1 %s% graph1 %u% (Z0-graph1)
  
  list(graph1=graph1,graph2=igraph::permute(graph2,permutation))
}

#' @export
#' @rdname sample_gnp
#' @return \code{sample_correlated_gnp_pair_w_junk} returns a list of two igraph object, named
#' \code{graph1} and \code{graph2}, which are two graphs whose adjacency matrix entries
#' correlated with \code{rho} and with first ncore vertices being core vertices and the rest being
#' junk vertices.
#' @examples
#' sample_correlated_gnp_pair_w_junk(50, 0.3, 0.5, 40)
#'
#'
sample_correlated_gnp_pair_w_junk <- function(n, corr, p, ncore=n,permutation=1:n,...){
  core <- 1:ncore
  junk <- (ncore+1):n
  
  cgnp_pair <- sample_correlated_gnp_pair(ncore,corr,p,...)
  
  if(ncore != n){
    pref_junk <- matrix(c(0,p,p,p),2,2)
    A <- sample_sbm(n,pref_junk,c(ncore,n-ncore),...)
    B <- sample_sbm(n,pref_junk,c(ncore,n-ncore),...)
    
    cgnp_pair <- with(cgnp_pair,{
      list(graph1=A %u% graph1,graph2=igraph::permute(B%u%graph2,permutation))
    })
  }
  cgnp_pair
}

