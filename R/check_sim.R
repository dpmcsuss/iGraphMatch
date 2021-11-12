#' @title Check the similarity matrix passed to a matching function
#'
#' @description Internal function that checks that a similarity matrix satisfies
#'  necessary conditions and modifies it for use in graph matching.
#'
#' @param sim Similarity matrix
#' @param seeds dataframe of seed matches from running \link{check_seeds}
#' @param nonseeds  dataframe of nonseed nodes from running \link{check_seeds}
#' @param totv1 total number of vertices in the first graph
#' @param totv2 total number of vertices in the second graph
#' @param for_nonseeds Whether the similarities are between non-seed nodes only (default = TRUE), or 
#'  if similarities among seed nodes are included (FALSE)
#'
#' @details The goal here is to be flexible in terms of the dimensions of the similarity matrix
#'  passed to \link{gm}. This is useful when the graphs have different orders in which case
#'  the function accepts matrices with dimensions equal to that of orders of the original graphs
#'  or the number of nonseeds.
#'
#' @return Standardized similarity matrix for similarities only between nonseeds across
#'  the two graphs, if for_nonseeds = TRUE, or between all nodes, if for_nonseeds = FALSE
#'
#' @rdname check_sim
check_sim <- function(sim, seeds, nonseeds, totv1, totv2, for_nonseeds = TRUE, square = TRUE){

  ns <- nrow(seeds)
  nn <- sapply(nonseeds, length)
  nv <- ns + nn

  # nv == max(totv1, totv2)

  # if its null then start with the zero matrix
  if (is.null(sim)) {
    sim <- Matrix::Matrix(0, nn[1], nn[2])
  }

  # if not we need to check dimensions
  dim_sim <- dim(sim)

  # TODO: it feels like this could get a little buggy based
  # on the fact the error checking is now a bit conditional on
  # square

  # first, if the sim is not square, we pad it to be square
  if(square && dim_sim[1] != dim_sim[2]){
    # has to be one of these dimensions
    if(all(dim_sim == c(totv1, totv2)) ||
        all(dim_sim + ns == c(totv1, totv2))){
      diff <- totv1 - totv2
      sim <- pad(sim, max(-diff, 0), max(diff, 0))
    } else {
      stop(paste0("Non square similarity matrices must have dimension equal to ",
        "that of the original graphs, ", totv1, " x ", totv2,
        ", or that of the nonseeds, ", totv1 - ns, " x ", totv2 - ns,
        "."))
    }

  }

  # now we've made them square if needed
  dim_sim <- dim(sim)

  # if they are nonseeds x nonseeds we're good
  if(for_nonseeds){
    if(dim_sim == nn || (square && dim_sim[1] == max(nn))){
      return(sim)
    } else {
      # otherwise keep only nonseeds
      return(sim[nonseeds$A, nonseeds$B])
    }
  } else {
    if(all(dim_sim + ns == c(totv1, totv2))) {
      return(pad(sim, ns, ns))
    }
    if(any(dim_sim < nv)){
      stop(paste0("Similarity matrices must have dimension equal to ",
                  totv1, " x ", totv2, "."))
    } else{
      return(sim)
    }
  }


  # otherwise, things seem wrong
  stop(paste0("Square similarity matrices must have dimension equal to the number of nonseeds, ",
      nn, ", or the total number of vertices, ", nv, "."))

}
