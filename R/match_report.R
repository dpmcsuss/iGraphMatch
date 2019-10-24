#' @title Matching performance summary
#'
#' @description Get a summary of the matching result and measures of the matching performance
#' based on several evaluation metrics associated with nodes and edges of two graphs.
#'
#' @param object A list. calls the matching result of applying a specific graph matching
#'   algorithm.
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param label A vector. NULL if the true correspondence between two graphs is unknown. 
#'   A vector indicating the true correspondence in the second graph if the true correspondence 
#'   is known, 
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
#' res <- graph_match_percolation(A, B, 1:4)
#' match_report(res, A, B)
#'
#' @export
#'
match_report <- function(object, A = A, B = B, label = NULL, ...){
  A <- A[]
  B <- B[]
  
  z <- object
  cat("Call: \n")
  print(z$call)
  
  # Matched nodes
  corr <- z$corr
  z$n.match <- nrow(corr) - z$ns
  cat("\n# Matches:", z$n.match)
  if(!is.null(label)){
    z$n.true.match <- sum(label[corr$corr_A] == corr$corr_B) - z$ns
    cat("\n# True Matches: ", z$n.true.match)
  }
  
  A_m <- A[corr$corr_A, corr$corr_A]
  B_m <- B[corr$corr_B, corr$corr_B]
  # Matched edges
  z$CE <- sum(A_m==B_m & A_m==1)
  z$CNE <- sum(A_m==B_m & A_m==0)
  z$EC <- z$CE/sum(A[]==1)
  cat("\n# Common Edges: ", z$CE,
      "\n# Common Non-edges: ", z$CNE)
  cat("\nEdge Correctness: ", z$EC)
  
  # objective value: ||A-PBP^T||_F
  z$Permutation <- get_perm(nrow(A), nrow(B), corr)
  z$Obj.Value <- Matrix::norm(A_m-B_m, type = "F")
  cat("\nObjective Value: ", z$Obj.Value)
  cat("\n")
}
