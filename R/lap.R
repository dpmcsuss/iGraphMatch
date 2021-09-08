

#' @title Linear (sum) assignment problem
#'
#' @description Compute the best bipartite matching
#' using one of three methods. For an n x n score matrix it find
#' \eqn{\max_{v\in \Pi_n} \sum_{i=1}^n score_{i, v(i)}}
#' where \eqn{\Pi_n} denotes all permutations on n objects.
#'
#' @param score matrix of pairwise scores
#' @param method One of "lapjv", "lapmod", "rect", or "clue"

#' @rdname do_lap
#'
#' @return \code{do_lap} returns a vector which indicates the
#'  best matching column for each row.
#'
#'
#' @details Solves a linear assignment using one of three methods.
#'  clue uses solve_lsap from the clue package.
#'  lapjv uses the Jonker-Volgenaut approach implemented in this package.
#'  lapmod use a version that exploits sparsity in the score matrix.
#'
#'
#' @examples
#' set.seed(12345)
#' cost <- Matrix::rsparsematrix(30, 30, .5)
#' cbind(
#'  do_lap(cost, "lapjv"),
#'  do_lap(cost, "lapmod"),
#'  do_lap(cost, "clue"),
#' )
#' do_lap(cost[1:5, ], "rect")
#' 
#' @export
do_lap <- function(score, method){
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
    # sinkhorn = {
    #   lambda <- 10
    #   n_iter <- 20
    #   sinkhorn(exp(lambda * score), n_iter)
    # },
    stop(paste0("The LAP method ", method,
        " is not implemented. Please use one of lapjv, lapmod, or clue."))
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


project_to_ds <- function(m, max_iter, tol) {
  # G_t <- [G_t + 1/n * (I - G_t + 11^T G_t/n)11^T - 11^T G_t/n)^+
}


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
  if(maximize) {
    x <- (1 + max(x)) - x
  }

  n <- nrow(x)
  m <- ncol(x)
  k <- n + m

  xpad <- matrix(max(x) + 1, k, k)
  xpad[1:n, 1:m] <- as.matrix(x)
  xpad[(n + 1):k, (m + 1):k] <- 0

  ind <- cpp_lapjv(xpad, maximize = FALSE)

  ind <- ind[1:n]
  ind[ind > m] <- -1

  return(ind)
}