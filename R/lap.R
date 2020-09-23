

#' @title Linear (sum) assignment problem
#'
#' @description Compute the best bipartite matching
#' using one of three methods. For an n x n score matrix it find
#' \eqn{\max_{v\in \Pi_n} \sum_{i=1}^n score_{i, v(i)}}
#' where \eqn{\Pi_n} denotes all permutations on n objects.
#'
#' @param score matrix of pairwise scores
#' @param method One off "lapjv", "lapmod", or "clue"
#' 
#'
#' @rdname do_lap
#'
#' @return \code{do_lap} returns a vector which indicates the 
#'  best matching column for each row.
#'  
#'
#' @details 
#' 
#' @examples
#' 
#' 
#' @export
do_lap <- function(score, method){
  n <- nrow(score)
  switch(method,
    lapjv = { 
      score <- as.matrix(score)
      rlapjv::lapjv(score, # round(score * n ^ 2 * max(score)),
        maximize = TRUE)
    },
    lapmod = {
      if( class(score) == "splrMatrix" ){
        rlapjv::lapmod(splr_to_sparse(score),
          maximize = TRUE)
      } else {
        rlapjv::lapmod(score, maximize = TRUE)
      }
    },
    clue = {
      score <- as.matrix(score)
      score <- score - min(score)
      as.vector(clue::solve_LSAP(score,
        maximum = TRUE))
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
  methods <- c("lapmod", "lapjv", "clue") #, "sinkhorn")
  if (!is.null(lap_method) && !(lap_method %in% methods)){
    stop(paste("Unrecognized LAP method:", lap_method,
      "Please use one of:", methods))
  }
  if (is.null(lap_method)){
    if("rlapjv" %in% rownames(utils::installed.packages()) &&
        totv1 / totv2 < 0.5 || totv2 / totv1 < 0.5){
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
