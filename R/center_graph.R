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
#' @param use_splr A boolean indicating whether to use the splrMatrix object when storing the 
#' centered graph.  Defaults to TRUE.
#'
#' @import irlba
#' 
#' @rdname center_graph
#' @return \code{center_graph} returns a centered adjacency matrix. 'Naive' scheme returns the
#' original adjacency matrix. 'Center' scheme returns a centered adjacency matrix with entries
#' equal to -1 or 1 where 1 corresponds to an edge. Scheme specified by a number returns a centered
#' adjacency matrix, which is calcalated by: A - Q, where Q is an approximate
#' estimation of the expectation of A by using Singular Value Thresholding.
#' @examples
#' A <- sample_correlated_gnp_pair(n = 10, corr = .5, p = .5)$graph1
#' center_graph(A, scheme = "naive")
#' center_graph(A, scheme = "center")
#' center_graph(A, scheme = 2)
#' center_graph(A, scheme = 1)
#'
#' @export
#'
center_graph <- function(A, scheme, use_splr = TRUE){
  if ( scheme == "naive" ){
    g <- A[]
  } else if ( scheme == "center" ){
    if (use_splr) {
      x <- 2*A[]
      g <- splr(x = x, a = rep(-1,dim(A[])[1]),b = rep(1,dim(A[])[2]))
    } else {
      g <- 2*A[] - 1
    }
  
  } else {
    r <- as.numeric(scheme)
    if (use_splr) {
      g <- splr(x=A[],a = -A[], rank = r, factorize = TRUE)
    } else {
      g <- A[] - low_rank_approx(A[], r)
    }
    
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



#' Pad a matrix object with extra rows/columns of 0s.
#' 
#' Attempts are made to make this padding efficient 
#' by employing sparse graphs
#' 
#' @param m matrix
#' @param nr number of rows to add
#' @param nc number of columns to add. (default = nr)
#' 
#' @returns m padded with nr rows and nc columns of zeros.
#' 
#' @export
pad <- function(m, nr, nc = nr){
  if(is(m, "splrMatrix")){
    dx <- as.integer(dim(m) + c(nr, nc))
    da <- dim(m@a)[2]
    splr(
      x = bdiag(m@x, Matrix(0, nr, nc)),
      a = Matrix(rbind2(m@a, Matrix(0, nr, da)), sparse = TRUE),
      b = Matrix(rbind2(m@b, Matrix(0, nc, da)), sparse = TRUE))
  }
  else{
    Matrix::bdiag(m, Matrix(0, nr, nc))
  }
}
