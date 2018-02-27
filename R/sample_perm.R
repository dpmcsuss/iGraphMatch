#' @title Sample random permutation matrix
#'
#' @description Sample an n-by-n random permutation matrix.
#'
#' @param n An integer. Dimension of the permutation matrix.
#'
#' @rdname sample_perm
#' @return \code{rperm} returns an n-by-n permutation matrix.
#' @examples
#' rperm(3)
#' @export
#'
rperm <- function(n){
  Matrix::Diagonal(n)[sample(n), ]
}
