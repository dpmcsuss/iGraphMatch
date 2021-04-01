#' @title Standardize seeds input data type
#'
#' @description Convert the input seeds data into data frame type with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}
#'
#' @param seeds A vector of integers or logicals, a matrix or a data frame. Input in the form of a
#' vector of integers denotes the indices of seeds which are identical in both graphs. Input in the
#' form of a vector of logicals indicate the location of seeds with TRUE and the indices of seeds
#' are identical in both graphs. Input in the form of a matrix or a data frame, with the first
#' column being the indices of \eqn{G_1} and the second
#' column being the corresponding indices of \eqn{G_2}.
#'
#' @param nv An integer. Number of total vertices.
#' @param logical An logical. TRUE indicates returns seeds in a vector of logicals where TRUE
#' indicates the corresponding vertex is a seed. FALSE indicates returns a data frame.
#'
#' @return returns a data frame with the first column being the corresponding indices of
#' \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2} or a vector of
#' logicals where TRUE indicates the corresponding vertex is a seed.
#'
#' @examples
#' #input is a vector of logicals
#' check_seeds(1:10 <= 3, nv = 10)
#'
#' #input is a vector of integers
#' check_seeds(c(1,4,2,7,3), nv = 10)
#'
#' #input is a matrix
#' check_seeds(matrix(1:4,2), nv = 10)
#'
#' #input is a data frame
#' check_seeds(as.data.frame(matrix(1:4,2)), nv = 10)
#'
#' @export
check_seeds <- function(seeds, nv, logical = FALSE){
  if(is.null(seeds)){
    seed_g1 <- numeric()
    seed_g2 <- numeric()
  }
  else if(is.logical(seeds)){
    seeds <- which(seeds==TRUE)
    seed_g1 <- seeds
    seed_g2 <- seeds
  } else if(is.vector(seeds)){
    seed_g1 <- seeds
    seed_g2 <- seeds
  } else if(is.matrix(seeds)){
    seed_g1 <- seeds[,1]
    seed_g2 <- seeds[,2]
  } else if(is.data.frame(seeds)){
    seeds <- as.matrix(seeds)
    seed_g1 <- seeds[,1]
    seed_g2 <- seeds[,2]
  } else{
    stop("Unrecognized seeds input format: seeds must be a vector of integers or logicals, a matrix or a data frame.")
  }

  if(logical==TRUE){
    seeds <- rep(FALSE, nv)
    seeds[seed_g1] <- TRUE
    seeds
  } else{
    seeds <- data.frame(A=seed_g1, B=seed_g2) # CORRECT this
    nonseeds <- data.frame(
        A = (1:nv)[!(1:nv %in% seeds$A)],
        B = (1:nv)[!(1:nv %in% seeds$B)])

    list(seeds = seeds, nonseeds = nonseeds)
  }

}
