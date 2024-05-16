#' @importFrom igraph sample_gnp sample_sbm graph_from_adjacency_matrix %u% %s% is_igraph V

#' @title Sample correlated G(n,p) random graphs
#'
#' @description Sample a pair of correlated G(n,p) random graphs with correlation between
#' two graphs being \code{corr} and edge probability being \code{p}.
#'
#' @param n An integer. Number of total vertices for the sampled graphs.
#' @param corr A number. The target Pearson correlation between the adjacency matrices
#' of the generated graphs. It must be in  [0,1] interval.
#' @param p A number. Edge probability between two vertices. It must be in open
#' [0,1] interval.
#' @param ncore An integer. Number of core vertices.
#' @param permutation A numeric vector to permute second graph.
#' @param ... Passed to \code{sample_gnp}.
#'
#' @rdname sample_gnp
#' @return \code{sample_correlated_gnp_pair} returns a list of two igraph object, named
#' \code{graph1} and \code{graph2}, whose adjacency matrix entries
#' are correlated with \code{corr}. If sample two graphs with junk vertices, the first
#' \code{ncore} vertices are core vertices and the rest are junk vertices.
#'
#' @references V. Lyzinski and D. E. Fishkind and C. E. Priebe (2014), \emph{Seeded Graph Matching
#' for Correlated Erdos-Renyi Graphs}.J. Mach. Learn. Res., pages 3513-3540.
#'
#' @examples
#' sample_correlated_gnp_pair(n=50, corr=0.3, p=0.5, ncore=40)
#' sample_correlated_gnp_pair(n=5, corr=0.3, p=0.5, permutation=c(1,3,2,4,5))
#'
#' @seealso \code{\link{sample_correlated_sbm_pair}}, \code{\link{sample_correlated_rdpg_pair}}
#'
#' @export
#'
sample_correlated_gnp_pair <- function(n, corr, p, ncore = n, permutation = 1:n, ...){
  if (p < 0 || p > 1) {
    stop("p must be between 0 and 1")
  }
  if (corr < 0 || corr > 1) {
    stop("corr must be between 0 and 1")
  }
  if (n <= 0) {
    stop("n must be positive integer")
  }
  if(ncore == n){
    sample_correlated_gnp_pair_no_junk(n, corr, p, permutation, ...)
  } else if(ncore < n){
    sample_correlated_gnp_pair_w_junk(n, corr, p, ncore, permutation, ...)
  } else{
    stop("ncore must be at most n.")
  }
}


sample_correlated_gnp_pair_no_junk <- function(n, corr, p, permutation=1:n, ...){

  # Make the first graph
  graph1 <- sample_gnp(n,p,...)

  # Make two graphs which will be used to make the
  # second graph
  Z0 <- sample_gnp(n,p*(1-corr),...)
  Z1 <- sample_gnp(n,p+corr*(1-p),...)

  graph2 <- Z1 %s% graph1 %u% (Z0-graph1)

  list(graph1=graph1,graph2=igraph::permute(graph2,permutation))
}


sample_correlated_gnp_pair_w_junk <- function(n, corr, p, ncore=n, permutation=1:n,...){
  core <- 1:ncore
  junk <- (ncore+1):n

  cgnp_pair <- sample_correlated_gnp_pair_no_junk(ncore,corr,p,...)

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

