#' @title Initialization of the start matrix
#'
#' @description Initialize the start matrix for graph matching iteration.
#'
#' @param start A matrix or a character. Any \code{nns-by-nns} doubly stochastic matrix or start
#' method like "bari", "convex" or "rds" to initialize the start matrix.
#' @param nns An integer. Number of non-seeds.
#' @param ns An integer. Number of seeds.
#' @param soft_seeds A vector, a matrix or a data frame. If there is no error in the soft seeds,
#' input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft seeds,
#' input should be in the form of a matrix or a data frame, with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}.
#' Note that if there are seeds in graphs, seeds should be put before non-seeds.
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}. Needed only when start is
#' convex.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}. Needed only when start is
#' convex.
#' @param seeds A logical vector. \code{TRUE} indicates the corresponding vertex is a seed. Needed
#' only when start is convex.
#' @param g A number. Specified in the range of [0, 1] to set weights to random permutaion matrix and
#' barycenter matrix.
#'
#' @rdname init_start
#' @return \code{init_start} returns a \code{nns-by-nns} doubly stochastic matrix as the start
#' matrix in the graph matching iteration. If conduct a soft seeding graph matching, returns a
#' \code{nns-by-nns} doubly stochastic matrix with 1's corresponding to the soft seeds and values
#' at the other places are derived by different start method.
#'
#' @examples
#' ss <- matrix(c(5, 4, 4, 3), nrow = 2)
#' # initialize start matrix without soft seeds
#' init_start(start = "bari", nns = 5)
#' init_start(start = "rds", nns = 3)
#' init_start(start = "rds_perm_bari", nns = 5)
#'
#' # initialize start matrix with soft seeds
#' init_start(start = "bari", nns = 5, ns = 3, soft_seeds = c(5, 7, 8))
#' init_start(start = "rds", nns = 5, soft_seeds = ss)
#' init_start(start = "rds_perm_bari", nns = 5, soft_seeds = ss)
#'
#' # initialize start matrix for convex graph matching
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:10 <= 2
#' init_start(start = "convex", nns = 8, A = g1, B = g2, seeds = seeds)
#'
#' # FW graph matching with incorrect seeds to start at convex start
#' init_start(start = "convex", nns = 8, ns = 2, soft_seeds = ss, A = g1, B = g2, seeds = seeds)
#'
#' # 
#'
#' @export
init_start <- function(start, nns, ns = 0, soft_seeds = NULL, ...){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  if (inherits(start, c("matrix", "Matrix"))){
    if (nss > 0 && any(dim(start) != c(nns - nss, nns - nss))) {
      stop("You are trying to use soft seeds but you've already
        specified the values for those rows and columns.")
    }
  } else if (is.function(start)) {
    sf <- start
    tryCatch(
      start <- sf(nns, ns, soft_seeds, ...),
      error = function(e){
        print(e)
        stop("Functions passed to init_start must have",
          "at least the arguments nns, ns, and softs_seeds")
      })

    # if we get back a full size matrix then just return
    if (all(dim(start) == nns)) {
      return(start)
    }
    # otherwise add in seeds below
  } else if (start == "bari"){
    start <- bari_start(nns - nss, ns)
  } else if (start == "rds") {
    start <- rds_sinkhorn_start(nns - nss, ns)
  } else if (start == "rds_perm_bari") {
    start <- rds_perm_bari_start(nns - nss, ns, ...)
  } else if (start == "convex") {
    # start at bari with soft seeds
    start <- init_start("bari", nns, ns, soft_seeds)
    # match and pull out doubly stochastic
    start <- graph_match_convex(..., start = start)$D 

    # don't add back in soft seeds b/c we've used them for convex
    # maybe message/warning about this being a silly thing
    # to do
    if (exists("seeds")) {
      nonseeds <- check_seeds(seeds, nv = nns + ns)$nonseeds
      start <- start[nonseeds$A, nonseeds$B]
    }
    return(start)
  } else {
    stop('Start must be either a matrix, function, or one of "bari", "rds",
      "rds_perm_bari", "convex"')
  }
  add_soft_seeds(start, nns, ns, soft_seeds)
}


add_soft_seeds <- function(start, nns, ns = 0, soft_seeds = NULL) {
  if (is.null(soft_seeds))
    return(start)

  cs <- check_seeds(soft_seeds, nv = nns + ns)
  seed_g1 <- cs$seeds$A - ns
  seed_g2 <- cs$seeds$B - ns
  cs <- check_seeds(cs$seeds - ns, nv = nns)
  nonseed_g1 <- cs$nonseeds$A
  nonseed_g2 <- cs$nonseeds$B

  new_start <- Matrix::Matrix(0, nns, nns)
  new_start[cbind(seed_g1, seed_g2)] <- 1
  new_start[nonseed_g1, nonseed_g2] <- start
  new_start
}

