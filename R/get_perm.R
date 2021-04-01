#' @title Get Permutation
#'
#' @description Get an \code{m-by-n} permutation matrix according to the mapping
#'   correspondence.
#'
#' @param match Either a graphMatch object or 2-column matrix or data frame.
#'  The first and second columns correspond to indices in \eqn{G_1} and 
#'  \eqn{G_2} respectively.
#' @param dim desired dimensions of the matrix. Note, this does not
#'  have to be square. If NULL and match is a graphMatch object then
#'  dim is set to dim(match)
#' @param padded If FALSE then this returns a square matrix the size of
#'  the larger of the two graph otherwise dim = dim(match). This is 
#'  ignored if match is not a graphMatch object.
#' @param seeds Whether to keep the seed vertices (TRUE) from the match
#'  or to remove them (FALSE). Ignored if match is not a graphMatch object.
#'  
#' 
#'
#' @rdname get_perm_mat
#'
#' @return \code{get_perm_mat} returns an \code{m-by-n} sparse permutation matrix or whose
#'   submatrix is a permutation matrix if only parts of nodes from both graphs get
#'   matched or in the case of matching graphs of different order.
#'
#' @examples
#' # returns a permutation matrix: m=n, all the nodes get matched
#' corr <- data.frame(corr_A = c(1,2,3,4), corr_B = c(1,4,2,3))
#' get_perm_mat(corr, c(4, 4))
#'
#' # submatrix is a permutation matrix: parts of graphs get matched
#' get_perm_mat(corr, c(5, 6))
#'
#' @export
#'
get_perm_mat <- function(match, dim = NULL, padded = FALSE, seeds = TRUE) {
  if (inherits(match, "graphMatch")) {
    if (is.null(dim) & padded) {
      # if we want padded then the perm should be as
      # large as the largest graph
      dim <- rep(max(length(match), dim(match)), 2)
      corr <- match@corr
    } else if(is.null(dim)){
      # otherwise the perm will be possibly rectangular
      # based on sizes of the two graphs
      dim <- dim(match)
   
      nmin <- min(dim(match), length(match))
      if (nmin == dim(match)[1]) {
        corr <- match[match[,1] %in% 1:nmin, ]
      } else if (nmin == dim(match)[1]) {
        corr <- match[match[,2] %in% 1:nmin, ]
      } else {
        corr <- match@corr
      }
    } else {
      corr <- match@corr
    }
  } else  {
    corr <- match
    if (is.null(dim)) {
      dim <- rep(0,2)
      dim[1] <- max(corr[,1])
      dim[2] <- max(corr[,2])
    } 
  }


  if (!is.null(dim)) {
    corr <- corr[corr[,1] <= dim[1] & corr[,2] <= dim[2], ]
    if (nrow(corr) == 0) {
      stop("The selected dim for the perm mat resulted in 0 valid matches being present")
    }
  }
  corr <- as.matrix(corr)
  perm <- Matrix::Matrix(0, dim[1], dim[2])
  perm[corr] <- 1

  if (!seeds & inherits(match, "graphMatch")) {
    nonseeds <- check_seeds(match$seeds, max(dim(match)))$nonseeds
    ns_a <- nonseeds[nonseeds[,1] <= dim[1], 1]
    ns_b <- nonseeds[nonseeds[,2] <= dim[2], 2]
    perm <- perm[ns_a, ns_b]
  }

  perm
}