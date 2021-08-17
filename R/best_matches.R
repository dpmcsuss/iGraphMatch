#' @title Find best matches
#'
#' @description Find a set of vertex-pairs in  order of a goodness of matching metric
#'
#' @param A A matrix, an \code{igraph} object, or a list of either. See \link{check_graph}
#' @param B A matrix, an \code{igraph} object, or a list of either. See \link{check_graph}
#' @param match \link{graphMatch}, eg result of call to \link{gm}
#' @param measure One of "row_cor", "row_diff", or "row_perm_stat". Measure for computing
#' goodness of matching.
#' @param num A positive integer. Number of pairs of best matched vertices needed.
#'
#' @return \code{best_matches} returns a data frame with the indices of best matched vertices
#' in \eqn{G_1} named \code{A_best}, the indices of best matched vertices in \eqn{G_2} named
#' \code{B_best} and the values of measure for best matches.
#'
#' \code{row_cor} takes 1 minus the row correlation value for the corresponding vertex.
#' \code{row_diff} takes the row difference value for each corresponding vertex.
#' \code{row_perm_stat} uses the row permutation statistics value.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 50, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:50 <= 10
#' nonseeds <- !seeds
#' match <- gm(g1, g2, seeds, method = "indefinite")
#'
#' # Application: select best matched seeds from non seeds as new seeds, and do the
#' # graph matching iteratively to get higher matching accuracy
#' best_matches(A = g1, B = g2, match = match, measure = "row_perm_stat", num = 5)
#'
#'
#' @export
#'
best_matches <- function(A, B, match, measure, num){

  if (!(measure %in% c("row_cor", "row_diff", "row_perm_stat"))) {
    stop('measure must be one of "row_cor", "row_diff", or "row_perm_stat".')
  }
  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  nv <- nrow(A[[1]])
  
  if (num < 0 || num > nv) {
    stop('num must be > 0 and <= number of nodes ', nv)
  }
  
  nc <- length(A)
  x <- !check_seeds(match$seeds, nv, logical = TRUE)
  match_corr <- match@corr
  x <- x[match_corr[,1]]
  match_corr <- match_corr[x,]
  stat <- rep(0, sum(x))

  for (ch in 1:nc) {
    A[[ch]] <- A[[ch]][match_corr[,1], match_corr[,1]]
    B[[ch]] <- B[[ch]][match_corr[,2], match_corr[,2]]

    # calculate measure stat
    stat <- stat + do.call(measure,list(A[[ch]],B[[ch]]))
  }

  # find top ranking nodes pairs
  if(measure != "row_cor"){
    stat <- -stat
  }
  rperm <- sample(length(stat))
  stat <- stat[rperm]
  if(num <= nv){
    topindex <- order(stat, decreasing = TRUE)[1:num]
  } else{
    stop("num can't exceed the total number of vertices.")
  }

  if(measure != "row_cor"){
    measure_value <- -stat[topindex]
  } else{
    measure_value <- stat[topindex]
  }
  top_matches <- match_corr[rperm[topindex],]

  best_matches <- data.frame(A_best=top_matches$corr_A, B_best=top_matches$corr_B,
                             measure_value)
  best_matches
}
