#' @title Sample graphs from edge probability matrix and correlation matrix
#'
#' @description Sample a pair of graphs with specified edge probability and
#' correlation between each pair of vertices.
#'
#' @param n An integer. Number of total vertices for the sampled graphs.
#' @param p_mat An \code{n-by-n} matrix. Edge probability matrix, each entry
#' should be in the open (0,1) interval.
#' @param c_mat An \code{n-by-n} matrix. The target Pearson correlation matrix,
#' each entry should be in the open (0,1) interval.
#' @param directed A logical. \code{TRUE} if the sampled graphs are directed.
#' @param X A matrix. Dot products matrix, each entry must be in open (0,1)
#' interval.
#' @param rho A number. The target Pearson correlation between the adjacency
#' matrices of the generated graphs. It must be in open (0,1) interval.
#' @param directed Logical scalar, whether to generate directed graphs.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @param permutation A numeric vector,permute second graph.
#' @param nc An integer. Number of core vertices.
#' @param ... Passed to \code{sample_correlated_rdpg_pair}.
#'
#' @rdname sample_ieg
#' @return \code{sample_correlated_ieg_pair} returns two igraph objects named
#' \code{graph1} and \code{graph2}.
#' @examples
#' n <- 50
#' p_mat <- matrix(runif(n^2),n)
#' c_mat <- matrix(runif(n^2),n)
#' sample_correlated_ieg_pair(n,p_mat,c_mat)
#'
#' @export
sample_correlated_ieg_pair<- function(n,p_mat,c_mat,directed=FALSE,loops=FALSE,permutation=1:n){
  g1 <- matrix(stats::rbinom(n^2,1,p_mat),n)
  z0 <- matrix(stats::rbinom(n^2,1,p_mat*(1-c_mat)),n)
  z1 <- matrix(stats::rbinom(n^2,1,p_mat*(1-c_mat)+c_mat),n)
  g2 <- z1*g1+z0*(1-g1)
  
  if(directed){
    mode <- "directed"
  } else{
    g1[row(g1)>=col(g1)] <- 0
    g2[row(g2)>=col(g2)] <- 0
    mode <- "undirected"
  }
  list(graph1 = graph_from_adjacency_matrix(g1, 
         mode = mode, diag = loops),
       graph2 = igraph::permute(graph_from_adjacency_matrix(g2,
         mode = mode, diag = loops),permutation))
}

#' @rdname sample_ieg
#' @return \code{sample_correlated_rdpg} returns two igraph objects named
#' \code{graph1} and \code{graph2} that are sampled from random dot product
#' graphs model.
#' @examples
#' ## sample a pair of igraph objects from random dot product graphs model
#' with dimension 3 and scale 8
#' n <- 50
#' xdim <- 3
#' scale <- 8
#' X <- matrix(rgamma(n*(xdim+1),scale,1),n,xdim+1)
#' X <- X/rowSums(X)
#' X <- X[,1:xdim]
#' sample_correlated_rdpg(X,rho=0.5)
#'
#' @export
sample_correlated_rdpg <- function(X,rho,nc=nrow(X),...){
  p_mat <- X%*%t(X)
  n <- nrow(X)
  if(length(rho)==1){
    c_mat <- matrix(0,n,n)
    c_mat[1:nc,1:nc] <- rho
  }else{
    c_mat <- rho
  }
  sample_correlated_ieg_pair(n,p_mat,c_mat,...)
}
