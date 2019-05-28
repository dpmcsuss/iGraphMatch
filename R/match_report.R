#' @title Matching report
#'
#' @description Get a report of matching performance based on several evaluation metrics
#'  according to the mapping correspondence.
#'
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param corr A matrix or a data frame. Matching correspondence with the first
#'   and second columns correspond to indices in \eqn{G_1} and \eqn{G_2} respectively.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   there is no error in seeds input can be a vector of seed indices in
#'   \eqn{G_1}. Or if there exists error in seeds, input in the form of a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param label A logical. TRUE if the true correspondence is known.
#'
#' @rdname match_report
#'
#' @return \code{match_report} returns a list of matching performance evaluation metrics
#' including number of matches, true matches, common edges, common non-edges, edge correctness
#' which is the fraction of common edges over number of edges in the first graph, and the
#' objective value ||A-PBP^T||_F.
#'
#' @examples
#' graphs <- sample_correlated_gnp_pair(10, .5, .3)
#' A <- graphs$graph1
#' B <- graphs$graph2
#' corr <- graph_match_percolation(A, B, 1:4)$corr
#' match_report(A, B, corr, seeds = 1:4)
#'
#' @export
#'
match_report <- function(A, B, corr, seeds = NULL, label = TRUE){
  A <- A[]
  B <- B[]
  seeds <- check_seeds(seeds)
  ns <- nrow(seeds)
  corr <- as.matrix(corr)

  # # matched pairs: not including seeds
  cat("# Matches: ", nrow(corr) - ns)
  # # true matches: not including seeds
  if(label == TRUE){
    cat("\n# True Matches: ", sum(corr[,1]==corr[,2]) - ns)
  }
  # # common edges
  CE <- sum(A==B & A==1)
  CNE <- sum(A==B & A==0)
  cat("\n# Common Edges: ", CE,
      "\n# Common Non-edges: ", CNE)

  # edge correctness: common edges / |E_1|
  cat("\nEdge Correctness: ", CE/sum(A==1))

  # objective value: ||A-PBP^T||_F
  P <- get_perm(nrow(A), nrow(B), corr)
  PB <- Matrix::crossprod(P, B)
  PBPT <- Matrix::tcrossprod(PB, P)
  cat("\nObjective Value: ", Matrix::norm(A - PBPT, type = "F"))
}


