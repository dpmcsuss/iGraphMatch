#' @title Convert soft seeds input data type
#'
#' @description Convert the soft seeds data into data frame type with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}
#'
#' @param soft_seeds A vector, a matrix or a data frame. If there is no error in soft
#' seeds, input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft
#' seeds, input in the form of a matrix or a data frame, with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}. Note
#' that if there are seeds in graphs, seeds should be put before non-seeds.
#'
#' @return returns a data frame with the first column being the indices of \eqn{G_1} and the
#' second column being the corresponding indices of \eqn{G_2}.
#'
#' @examples
#' ##input is a vector
#' check_seeds(c(1,4,2,7,3))
#'
#' ##input is a matrix
#' check_seeds(matrix(1:4),2)
#'
#' ##input is a data frame
#' check_seeds(as.data.frame(matrix(1:4),2))
#'
#' @export
check_seeds <- function(soft_seeds){
  if(is.vector(soft_seeds)){
    seed_g1 <- soft_seeds
    seed_g2 <- soft_seeds
  } else if(is.matrix(soft_seeds)){
    seed_g1 <- soft_seeds[,1]
    seed_g2 <- soft_seeds[,2]
  } else{
    soft_seeds <- as.matrix(soft_seeds)
    seed_g1 <- soft_seeds[,1]
    seed_g2 <- soft_seeds[,2]
  }

  data_frame(seed_A=seed_g1, seed_B=seed_g2)
}
