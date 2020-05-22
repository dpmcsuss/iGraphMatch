
check_sim <- function(sim, seeds, nonseeds, totv1, totv2){

  ns <- nrow(seeds)
  nn <- nrow(nonseeds)
  nv <- ns + nn

  # nv == max(totv1, totv2)

  # if its null then return the zero matrix
  if(is.null(sim)){
    return(Matrix::Matrix(0, nn, nn))
  }

  # if not we need to check dimensions
  dim_sim <- dim(sim)

  # first, if the sim is not square, we pad it to be square
  if(dim_sim[1] != dim_sim[2]){
    # has to be one of these dimensions
    if( all(dim_sim == c(totv1, totv2)) || 
        all(dim_sim + ns == c(totv1, totv2)) ){
      diff <- totv1 - totv2
      sim <- pad(sim, max(-diff, 0), max(diff, 0))
    } else {
      stop(paste0("Non square similarity matrices must have dimension equal to ",
        "that of the original graphs, ", totv1, " x ", totv2,
        ", or that of the nonseeds, ", totv1 - ns, " x ", totv2 - ns,
        "."))
    }

  }

  # now we've made them square
  dim_sim <- dim(sim)[1]

  # if they are nonseeds x nonseeds we're good
  if(dim_sim == nn){
    return(sim)
  } else if(dim_sim == nv){
    # otherwise keep only nonseeds
    return(sim[nonseeds$A, nonseeds$B])
  } 

  # otherwise, things seem wrong
  stop(paste0("Similarity matrix must be either NULL or ", 
      "a square matrix of dimension equal to the number of nonseeds, ", 
      nn, ", or the total number of vertices, ", nv, "."))

}