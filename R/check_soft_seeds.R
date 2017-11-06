#' @title Check initialization for soft seeding
#'
#' @description Initialize the start matrix for soft seeding graph matching iteration.
#'
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or character value
#' like "bari", "convex" or "rds" to initialize the starting matrix.
#' @param nns An integer. Number of non-seeds.
#' @param ns An integer. Number of seeds.
#' @param soft_seeds A vector, a matrix or a data frame. If there is no error in soft
#' seeds, input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft
#' seeds, input in the form of a matrix or a data frame, with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}. Note
#' that if there are seeds in graphs, seeds should be put before non-seeds.
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}. Needed only when start is
#' convex.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}. Needed only when start is
#' convex.
#' @param seeds A logical vector. \code{TRUE} indicates the corresponding vertex is a seed. Needed
#' only when start is convex.
#' @param non_seed_core A logical vector. \code{TRUE} indicates the corresponding vertex is a
#' non_seed_core vertex. Needed only when start is convex.
#'
#' @rdname check_soft_seeds
#' @return \code{check_soft_seeds} returns a \code{nns-by-nns} matrix with 1's corresponding to the
#' adaptive seeds and being bari-centered at other places.
#'
#' @examples
#' ss <- matrix(c(2,4,2,3), nrow = 2)
#' check_soft_seeds("bari", 5, 3, c(5,7,8))
#' check_soft_seeds(start = "rds", nns = 5, soft_seeds = ss)
#'
#' cgnp_pair <- sample_correlated_gnp_pair(n = 50, rho =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:50 <= 10
#' ns_core <- 1:50 <= 50
#' check_soft_seeds(start = "convex", nns = 5, soft_seeds = ss, A = g1, B = g2, seeds = seeds, non_seed_core = ns_core)
#'
#' @export
check_soft_seeds <- function(start, nns, ns = 0, soft_seeds, A = NULL, B = NULL, seeds = NULL, non_seed_core = NULL){
  if(grepl("atrix",class(start))){
    start <- start
  } else if(start == "bari"){
    start <- bari_start(nns,ns,soft_seeds)
  } else if(start =="rds"){
    start <- rds_sinkhorn_start(nns,ns,soft_seeds)
  } else if(start == "convex"){
    D <- graph_match_adpt_seeds(A, B, seeds, non_seed_core, soft_seeds, start="convex")$D_hard
    start <- D[!seeds,!seeds]
  }
  start
}
