#' @title Sample graphs pair from stochastic block model
#'
#' @description Sample a pair of random graphs from stochastic block model with correlation between
#' two graphs being \code{rho} and edge probability being \code{p}.
#'
#' @param n An integer. Number of vertices in the graph.
#' @param pref.matrix The matrix giving the Bernoulli rates. This is a \code{K-by-K} matrix, where
#' \code{k} is the number of groups. The probability of creating an edge between vertices from groups
#' \code{i} and \code{j} is given by element \code{i,j}. For undirected graphs, this matrix must be
#' symmetric.
#' @param block.sizes A numeric vector. Give the number of vertices in each group. The sum of the
#' vector must match the number of vertices.
#' @param rho A number. The target Pearson correlation between the adjacency matrices of the generated 
#' graphs. It must be in open (0,1) interval.
#' @param permutation. A numeric vector, permute second graph. 
#' @param core.block.size A numeric vector. Give the number of core vertices in each group. Entries
#' should be smaller than \code{block.sizes} and the vector length should be the same as \code{block.sizes}.
#' @param ... Passed to \code{sample_correlated_sbm_pair}.
#'
#' @rdname sample_sbm
#' @return A list of two igraph object, named \code{graph1} and \code{graph2}.
#' @examples
#' pm <- cbind( c(.1, .001), c(.001, .05) )
#' sample_correlated_sbm_pair(1000, pref.matrix=pm, block.sizes=c(300,700), rho=0.5)
#' @export
#'
sample_correlated_sbm_pair <- function(n, pref.matrix, block.sizes, rho, permutation=1:n, ...){
  
  K <- length(block.sizes)
  # Make the first graph
  graph1 <- sample_sbm(n,pref.matrix,block.sizes,...)
  
  # Make two graphs which will be used to make the
  # second graph
  corr.matrix <- (1-rho)*pref.matrix
  Z0 <- sample_sbm(n,corr.matrix,block.sizes,...)
  Z1 <- sample_sbm(n,corr.matrix+rho,block.sizes,...)
  
  graph2 <- Z1 %s% graph1 %u% (Z0-graph1)
  
  list(graph1=graph1,graph2=permute(graph2,permutation))
}

#' @rdname sample_sbm
#' @examples
#' sample_correlated_sbm_pair_w_junk(1000, pref.matrix=pm, block.sizes=c(300,700), rho=0.5,
#' core.block.sizes=c(200,500))
#' @export
#'
sample_correlated_sbm_pair_w_junk <- function(
  n, pref.matrix, block.sizes, rho, core.block.sizes, permutation=1:n, ...){
  
  K <- length(block.sizes)
  ncore <- sum(core.block.sizes)
  core <- 1:ncore
  junk <- (ncore+1):n
  
  junk.block.sizes <- block.sizes - core.block.sizes
  all.block.sizes <- c(core.block.sizes,junk.block.sizes)
  all.pref.matrix <- kronecker(matrix(1,2,2),pref.matrix)
  # Make the first graph
  graph1 <- sample_sbm(n,all.pref.matrix,all.block.sizes,...)
  
  # Make two graphs which will be used to make the
  # second graph
  all.pref.matrix <- kronecker(matrix(c(1-rho,1,1,1),2),pref.matrix)
  Z0 <- sample_sbm(n,all.pref.matrix,all.block.sizes,...)
  
  all.pref.matrix[1:K,1:K] <- all.pref.matrix[1:K,1:K] + rho
  Z1 <- sample_sbm(n,all.pref.matrix,all.block.sizes,...)
  
  graph2 <- Z1 %s% graph1 %u% (Z0-graph1)
  
  list(graph1=graph1,graph2=permute(graph2,permutation))
}

