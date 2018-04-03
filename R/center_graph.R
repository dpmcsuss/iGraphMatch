#' @title Center adjacency matrix
#'
#' @description Center the adjacency matrix, including center the adjacency matrix to entries
#' equal to -1 or 1, center the adjacency matrix by using Universal Singular Value Thresholding.
#'
#' @param A A matrix or an igraph object. Adjacency matrix.
#' @param scheme A number or a character. Type of matrix to return. If approximatedly center
#' the adjacency matrix using Universal Singular Value Thresholding, specify \code{scheme} with
#' the number of singular values used to estimate the expectation of adjacency matrix. Suitable
#' choice of scheme can be obtained with rank of order logn, where n is the dimension of A, expecially
#' in the setting of latent space graph models.
#'
#' @rdname center_graph
#' @return \code{center_graph} returns a centered adjacency matrix. 'Naive' scheme returns the
#' original adjacency matrix. 'Center' scheme returns a centered adjacency matrix with entries
#' equal to -1 or 1 where 1 corresponds to an edge. Scheme specified by a number returns a centered
#' adjacency matrix, which is calcalated by: A - \hat{Q}, where \hat{Q} is an approximate
#' estimation of the expectation of A by using Universal Singular Value Thresholding.
#' @examples
#' A <- sample_correlated_gnp_pair(n = 10, corr = .5, p = .5)$graph1
#' center_graph(A, scheme = "naive")
#' center_graph(A, scheme = "center")
#' center_graph(A, scheme = 2)
#'
#' @export
#'
center_graph <- function(A, scheme){
  if ( scheme == "naive" ){
    g <- A[]
  } else if ( scheme == "center" ){
    g <- 2*A[] - 1
  } else {
    r <- as.numeric(scheme)
    g <- A[] - low_rank_approx(A[], r)
  }
  g
}

low_rank_approx <- function(A,ndim){
  usv <- irlba(A, ndim)
  if ( ndim > 1 ){
    with(usv, u %*% diag(d) %*% t(v))
  } else{
    with(usv, outer(as.vector(u), as.vector(v)) * d)
  }
}
