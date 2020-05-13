

check_graph <- function(A, B, same_order = TRUE, square = TRUE, as_list = TRUE, as_igraph = FALSE){

  # **********NEED to implement PARAMETERS**********
  # Right now defaults are only option

  if(as_igraph){
    stop("Not yet implemented")
    # return(check_graph_igraph(A, B, same_order))
  }

  # this will make the graphs be matrices if they are igraph objects
  if(is.list(A) && !igraph::is.igraph(A)){
    A <- lapply(A, function(Al) Al[])
  } else {
    A <- list(A[])
  }
  if( is.list(B) && !igraph::is.igraph(B)){
    B <- lapply(B, function(Bl) Bl[])
  } else {
    B <- list(B[])
  }

  totv1 <- ncol(A[[1]])
  totv2 <- ncol(B[[1]])

  if(any(sapply(A, function(Al) ncol(Al) != totv1))){
    stop("A contains graphs of different orders. For multiple graph matching, all graphs must have the same number of vertices.")
  }
  if(any(sapply(B, function(Bl) ncol(Bl) != totv2))){
    stop("B contains graphs of different orders. For multiple graph matching, all graphs must have the same number of vertices.")
  }
  # Check for square
  if(any(sapply(A, function(Al) nrow(Al) != totv1))){
    stop("A is not square. graph_match_FW only supports ",
      "square matrices for matching.")
  }
  if(any(sapply(B, function(Bl) nrow(Bl) != totv2))){
    stop("B is not square. graph_match_FW only supports ",
      "square matrices for matching.")
  }

  try({
    A <- lapply(A, function(Al) as(Al, "dgCMatrix"))
    B <- lapply(B, function(Bl) as(Bl, "dgCMatrix"))
  }, silent = TRUE)
  try({B <- as(B, "dgCMatrix")}, silent = TRUE)

  if(same_order){
    if(totv1 > totv2){
      diff <- totv1 - totv2
      B <- lapply(B, function(Bl)
        pad(Bl[], diff))
    }else if(totv1 < totv2){
      diff <- totv2 - totv1
      A <- lapply(A, function(Al)
        pad(Al[], diff))
    }
  }

  if(! as_list){
    if(length(A) > 1){
      warning("A is multi-layer and must be converted to single layer.")
    } else if(length(A) > 1) {
      warning("B is multi-layer and must be converted to single layer.")
    } else {
      A <- A[[1]]
      B <- B[[1]]
    }
  }

  list(g1 = A, g2 = B, totv1 = totv1, totv2 = totv2)
}