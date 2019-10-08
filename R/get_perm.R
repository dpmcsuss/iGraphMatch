#' @title Get Permutation
#'
#' @description Get an $m$-by-$n$ permutation matrix according to the mapping
#'   correspondence.
#'
#' @param m An integer. Order of \eqn{G_1}.
#' @param n An integer. Order of \eqn{G_2}.
#' @param corr A matrix or a data frame. Matching correspondence with the first
#'   and second columns correspond to indices in \eqn{G_1} and \eqn{G_2} respectively.
#'
#' @rdname get_perm
#'
#' @return \code{get_perm} returns an $m$-by-$n$ sparse permutation matrix or whose
#'   submatrix is a permutation matrix if only parts of nodes from both graphs get
#'   matched or in the case of matching graphs of different order.
#'
#' @examples
#' # returns a permutation matrix: m=n, all the nodes get matched
#' corr <- data.frame(corr_A = c(1,2,3,4), corr_B = c(1,4,2,3))
#' get_perm(4, 4, corr)
#'
#' # submatrix is a permutation matrix: parts of graphs get matched
#' get_perm(5, 6, corr)
#'
#' @export
#'
get_perm <- function(m, n, corr){
  corr <- as.matrix(corr)
  perm <- Matrix::Matrix(0, m, n)
  perm[corr] <- 1
  perm
}
