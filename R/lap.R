

#' @title Linear (sum) assignment problem
#'
#' @description Compute the best bipartite matching
#' using one of three methods. For an n x n score matrix it find
#' \eqn{\max_{v\in \Pi_n} \sum_{i=1}^n score_{i, v(i)}}
#' where \eqn{\Pi_n} denotes all permutations on n objects.
#'
#' @param score matrix of pairwise scores
#' @param method One of "lapjv", "lapmod", or "clue"
#'
#'
#' @rdname do_lap
#'
#' @return \code{do_lap} returns a vector which indicates the
#'  best matching column for each row.
#'
#'
#' @details Solves a linear assignment using one of three methods.
#'  "clue" uses \code{solve_lsap} from the clue package.
#'  "lapjv" uses the Jonker-Volgenaut approach implemented in this package.
#'  "lapmod" use a modification of JV that exploits sparsity in the 
#'  score matrix.
#'  "rect" uses lapjv on a reduced problem. If m is the order of
#'  the smaller graph and n is the order of the larger graph, "rect"
#'  can yield significant computational time improvements if 
#'  $m = o(\sqrt{n})$.
#'
#'  Scores do not need to be non-negative. For "clue" the scores are pre-translated to be
#'  non-negative which preserves the LAP solution.
#'
#'
#' @references R. Jonker, A. Volgenant (1987). \emph{A shortest augmenting path algorithm
#' for dense and sparse linear assignment problems}. Computing, pages 325-340.
#' @references A. Volgenant (1996). \emph{Linear and Semi-Assignment Problems: A
#'   Core Oriented Approach}. Computer Ops Res., pages 917-932.
#' @references C. H. Papadimitriou and K. Steiglitz (1998). \emph{Combinatorial Optimization:
#' Algorithms and Complexity}. Courier Corporation.
#' @examples
#' set.seed(12345)
#' cost <- Matrix::rsparsematrix(10, 10, .5)
#' cbind(
#'  do_lap(cost, "lapjv"),
#'  do_lap(cost, "lapmod"),
#'  do_lap(cost, "clue")
#' )
#'
#' @export
do_lap <- function(score, method = "clue"){
  if (!inherits(score, c("matrix", "Matrix"))) {
    stop("score must a matrix-like object.")
  }
  diff <- dim(score)[1] - dim(score)[2]
  if (method != "rect" && diff != 0){
    score <- pad(score, max(-diff, 0), max(diff, 0))
  }
  n <- nrow(score)
  method <- set_lap_method(method, n, n)
  switch(method,
    lapjv = {
      score <- as.matrix(score)
      lapjv(score, # round(score * n ^ 2 * max(score)),
        maximize = TRUE)
    },
    lapmod = {
      if( class(score) == "splrMatrix" ){
        lapmod(splr_to_sparse(score),
          maximize = TRUE)
      } else {
        lapmod(score, maximize = TRUE)
      }
    },
    clue = {
      score <- as.matrix(score)
      score <- score - min(score)
      as.vector(clue::solve_LSAP(score,
        maximum = TRUE))
    },
    rect = {
      rect_lap(score)
    },
    stop(paste0("The LAP method '", method,
        "' is not implemented. Please use one of 'lapjv', 'lapmod',",
        " 'clue', or 'rect'."))
  )
}

set_lap_method <- function(lap_method, totv1, totv2){
  methods <- c("lapmod", "lapjv", "clue", "rect") #, "sinkhorn")
  if (!is.null(lap_method) && !(lap_method %in% methods)){
    stop(paste("Unrecognized LAP method:", lap_method,
      "Please use one of:", paste(methods, collapse = " ")))
  }
  if (is.null(lap_method)){
    if (totv1 / totv2 < 0.5 || totv2 / totv1 < 0.5) {
      lap_method <- "lapmod"
    } else {
      lap_method <- "clue"
    }
  }
  lap_method
}


# project_to_ds_l2 <- function(m, max_iter, tol) {
#   # Iterate G_t <- [G_t + 1/n * (I - G_t + 11^T G_t/n)11^T - 11^T G_t/n)^+
# }



rect_lap <- function(g_rect) {
  n = nrow(g_rect)
  m = ncol(g_rect)

  topk       = t(sapply(1:n, function(i) order(-g_rect[i,])))
  topk       = topk[,1:n]
  sel        = sort(unique(c(topk)))
  g_rect_sub = g_rect[,sel]

  # solve sub-problem
  ind        = sel[rect_lap_inner(g_rect_sub, maximize=TRUE)]

  # pad ind to get full solution
  ind_extra = 1:m
  ind_extra = ind_extra[!(ind_extra %in% ind)]
  ind       = c(ind, ind_extra)
  return(ind)
}


rect_lap_inner <- function(x, maximize) {
  x <- as.matrix(x)
  if(maximize) {
    x <- (1 + max(x)) - as.matrix(x)
  }

  n <- nrow(x)
  m <- ncol(x)
  k <- max(n, m)

  xpad <- matrix(max(x) + 1, k, k)
  xpad[1:n, 1:m] <- as.matrix(x)
  # xpad[(n + 1):k, (m + 1):k] <- 0

  ind <- cpp_lapjv(xpad, maximize = FALSE)

  ind <- ind[1:n]
  ind[ind > m] <- -1

  return(ind)
}
