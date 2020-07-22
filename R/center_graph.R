#' @title Center adjacency matrix
#'
#' @description Center the adjacency matrix, including center the adjacency matrix to entries
#' equal to -1 or 1, center the adjacency matrix by using Universal Singular Value Thresholding.
#'
#' @param A A matrix or an igraph object. Adjacency matrix.
#' @param scheme A character vector, number or pair of numbers. Default c(-1, 1). See Details.
#' 
#' @param use_splr A boolean indicating whether to use the splrMatrix object when storing the 
#' centered graph.  Defaults to TRUE.
#'
#' @details  The options for scheme are
#' \itemize{
#'  \item "naive" Returns original A
#'  \item Integer: Returns \eqn{A - A_{scheme}} where
#'    \eqn{A_{scheme}} is the best rank-scheme approximation
#'    of A.
#'  \item A pair of scalars: Returns s * A + a such that the
#'    minimum of the returned matrix is min(scheme) and the
#'    maximum is max(scheme).
#'  \item "center": Same as scheme=c(-1,1)
#' }
#' 
#' 
#' @rdname center_graph
#' @return centered adjacency matrix as a splrMatrix if
#'  useSplr = TRUE, otherwise as a Matrix object.
#' @examples
#' A <- sample_correlated_gnp_pair(n = 10, corr = .5, p = .5)$graph1
#' center_graph(A, scheme = "naive")
#' center_graph(A, scheme = "center")
#' center_graph(A, scheme = 2)
#' center_graph(A, scheme = c(-4, 2))
#'
#' @export
#'
center_graph <- function(A, scheme = c(-1, 1), use_splr = TRUE){
  if (!(length(scheme) %in% c(1,2))) {
    stop("scheme must be either 'center', 'naive', ",
      "a positive integer, or a pair of scales.")
  }

  if (is.character(scheme) && scheme == "center") {
    scheme <- c(-1, 1)
  }
  if (is.character(scheme) && scheme == "naive" ){
    g <- A[]
  } else if (length(scheme) == 1 && is.numeric(scheme)) {
    r <- as.integer(scheme)
    if (use_splr) {
      g <- splr(x = A[], a = -A[], 
        rank = r, factorize = TRUE)
    } else {
      g <- A[] - low_rank_approx(A[], r)
    }
    
  } else if (is.numeric(scheme) && length(scheme) == 2) {
    nmx <- max(scheme)
    nmn <- min(scheme)
    mx <- max(A[])
    mn <- min(A[])
    s <- (nmx - nmn) / (mx - mn)
    a <- - (mx + mn) / (mx - mn)  * (nmx - nmn) / 2 + 
        (nmx + nmn) / 2
    if (use_splr) {
      g <- splr(x = s * A[], 
        a = rep(a, dim(A[])[1]),
        b = rep(1, dim(A[])[2]))
    } else {
      g <- s * A[] + a
    } 
  } else {
    stop("scheme must be either 'center', 'naive', ",
      "a positive integer, or a pair of scalars.")
  }
  g
}


low_rank_approx <- function(A,ndim){
  usv <- irlba::irlba(A, ndim)
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
  else if(is(m, "matrix_list")) {
    m <- lapply(m, function(ml) ml[])
    matrix_list(lapply(m, pad, nr = nr, nc = nc))
  }
  else{
    Matrix::bdiag(m, Matrix(0, nr, nc))
  }
}
