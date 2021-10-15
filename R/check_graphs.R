
#' @title Parameter checking for a graph-pair
#'
#' @description Internal function that checks that the pair of graphs passed
#'  to a matching-related functions satisfies necessary conditions and modifies
#'  them according to specified parameters. check_single_graph
#'  does similar checks and modifications but just for one graph or list of graphs.
#'
#' @param A A matrix, an igraph object, or list of either.
#' @param B A matrix, an igraph object, or list of either.
#' @param same_order Whether the returned objects should have the same number of nodes.
#'  If the graphs start with different numbers of nodes the smaller graph is padded with
#'  isolated vertices. (default = TRUE)
#' @param square Whether the matrices need to be square. (default = TRUE)
#'  Currently non-square matrices are not supported.
#' @param as_list Whether to return the results as a matrix_list. (default = TRUE)
#'  If FALSE and A and B have length > 1
#' @param as_igraph Whether to return an igraph object. (default=FALSE)
#'  Only allowed if the original parameters are igraph objects.
#'  If FALSE, then this converts the objects to sparse matrices.
#'
#' @details If A and B are lists of matrices or igraph objects, then the lists
#'  must be the same length. Additionally, within each list the graphs need to have the
#'  same number of vertices but this does not need to be true across lists.
#'
#' @rdname check_graph
#' @return List containing A and B modified according to the parameters and the number of 
#'  vertices in each graph in totv1 and totv2.
#'  
#'
#' @export
#'
check_graph <- function(A, B,
  same_order = TRUE, square = TRUE, 
  as_list = TRUE, as_igraph = FALSE) {


  if (as_igraph) {
    if (igraph::is.igraph(A) && igraph::is.igraph(B)) {
        totv1 <- igraph::gorder(A)
        totv2 <- igraph::gorder(B)
        if (totv1 > totv2 && same_order) {
          B <- igraph::add_vertices(B, totv1 - totv2)
        }
        if (totv2 > totv1 && same_order) {
          A <- igraph::add_vertices(A, totv2 - totv1)
        }
        return(list(g1 = A, g2 = B,
          totv1 = totv1, totv2 = totv2))

    }
    stop("Check graph only supports as_igraph = TRUE if both A and B are igraph objects.")
    # return(check_graph_igraph(A, B, same_order))
  }

  # this will make the graphs be matrices if they are igraph objects
  if (is.list(A) && !igraph::is.igraph(A)) {
    A <- matrix_list(lapply(A, function(Al) Al[]))
  } else {
    A <- matrix_list(list(A[]))
  }
  if ( is.list(B) && !igraph::is.igraph(B)) {
    B <- matrix_list(lapply(B, function(Bl) Bl[]))
  } else {
    B <- matrix_list(list(B[]))
  }

  totv1 <- ncol(A[[1]])
  totv2 <- ncol(B[[1]])

  if (any(sapply(A, function(Al) ncol(Al) != totv1))) {
    stop("A contains graphs of different orders. For multiple graph matching, all graphs must have the same number of vertices.")
  }
  if (any(sapply(B, function(Bl) ncol(Bl) != totv2))) {
    stop("B contains graphs of different orders. For multiple graph matching, all graphs must have the same number of vertices.")
  }
  # Check for square
  if (square) {
    if (any(sapply(A, function(Al) nrow(Al) != totv1))) {
      stop("A is not square. This method only supports ",
        "square matrices for matching.")
    }
    if (any(sapply(B, function(Bl) nrow(Bl) != totv2))) {
      stop("B is not square. This method only supports ",
        "square matrices for matching.")
    }
  } else {
    stop("square = FALSE is not yet supported for check_graph")
  }


  try({
    A <- matrix_list(lapply(A, function(Al) as(Al, "dgCMatrix")))
    B <- matrix_list(lapply(B, function(Bl) as(Bl, "dgCMatrix")))
  }, silent = TRUE)
  # try({
  #   A <- as(A, "dgCMatrix")
  #   B <- as(B, "dgCMatrix")
  # }, silent = TRUE)

  if (same_order) {
    if (totv1 > totv2) {
      diff <- totv1 - totv2
      B <- pad(B, diff)
      # B <- lapply(B, function(Bl)
      #   pad(Bl[], diff))
    }else if (totv1 < totv2) {
      diff <- totv2 - totv1
      A <- pad(A, diff)
      # A <- lapply(A, function(Al)
      #   pad(Al[], diff))
    }
  }

  if (!as_list) {
    if (length(A) > 1) {
      stop("A is multi-layer and must be converted to single layer.\
       (check_graph: is_list = FALSE)")
    } else if (length(B) > 1) {
      stop("B is multi-layer and must be converted to single layer.\
       (check_graph: is_list = FALSE)")
    } else {
      A <- A[[1]]
      B <- B[[1]]
    }
  }

  list(g1 = A, g2 = B, totv1 = totv1, totv2 = totv2)
}


#' @rdname check_graph
check_single_graph <- function(A, square = TRUE, 
  as_list = TRUE, as_igraph = FALSE) {


  if (as_igraph) {
    if (igraph::is.igraph(A)) {
        return(A)
    }
    stop("Check single graph only supports as_igraph = TRUE if both A is an igraph object")
    # return(check_graph_igraph(A, B, same_order))
  }

  # this will make the graphs be matrices if they are igraph objects
  if (is.list(A) && !igraph::is.igraph(A)) {
    A <- matrix_list(lapply(A, function(Al) Al[]))
  } else {
    A <- matrix_list(list(A[]))
  }
  
  totv <- ncol(A[[1]])
  
  if (any(sapply(A, function(Al) ncol(Al) != totv))) {
    stop("A contains graphs of different orders. All layers must have the same number of vertices.")
  }
  if (square) {
    if (any(sapply(A, function(Al) nrow(Al) != totv))) {
      stop("A is not square. This method only supports ",
        "square matrices for matching.")
    }
  } else {
    stop("square = FALSE is not yet supported for check_graph")
  }


  try({
    A <- matrix_list(lapply(A, function(Al) as(Al, "dgCMatrix")))
  }, silent = TRUE)
  # try({
  #   A <- as(A, "dgCMatrix")
  #   B <- as(B, "dgCMatrix")
  # }, silent = TRUE)

  if (!as_list) {
    if (length(A) > 1) {
      stop("A is multi-layer and must be converted to single layer.\
       (check_graph: is_list = FALSE)")
    }  else {
      A <- A[[1]]
    }
  }

  A
}