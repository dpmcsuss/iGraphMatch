#' @title Initialization of the start matrix
#'
#' @description Initialize the start matrix for graph matching iteration.
#'
#' @param start A matrix, character, or function. A \code{nns-by-nns} matrix, start
#' method like "bari", "convex" or "rds", or a function  to initialize the start matrix.
#' If a function, it must have at least the arguments nns, ns, and softs_seeds.
#' @param nns An integer. Number of non-seeds.
#' @param ns An integer. Number of seeds.
#' @param soft_seeds A vector, a matrix or a data frame indicating entries of the start matrix
#'  that will be initialized at 1 to indicate . See \link{check_seeds}.
#' @param seeds  A vector, a matrix or a data frame. Indicating hard seeds.
#' These are used for "convex" start but otherwise are ignored. 
#' @param ... Arguments passed to other start functions. See details in Values section.
#'
#' @rdname init_start
#' @return \code{init_start} returns a \code{nns-by-nns} doubly stochastic matrix as the start
#' matrix in the graph matching iteration. If conduct a soft seeding graph matching, returns a
#' \code{nns-by-nns} doubly stochastic matrix with 1's corresponding to the soft seeds and values
#' at the other places are derived by different start method.
#'
#'@details
#' When \code{start} is a character, there are five options.
#' \itemize{
#' \item \code{"bari"} initializes at the barycenter.
#' \item \code{"rds_perm_bari"} gives a random linear combination of barycenter and
#'   a random permutation matrix, (1-a) B + a P. The argument \code{g} controls a
#'   with a being sampled as \code{g * runif()}.
#' \item \code{"rds"} gives a random doubly stochastic matrix. Users can specify a
#'   random deviates generator to the \code{distribution} argument, and the default is \code{runif}. 
#'   A random matrix with iid entries from \code{distribution} and the the Sainkhorn algorithm is applied
#'   to produce the output.
#' \item \code{"rds_from_sim"} gives a random doubly stochastic matrix derived from
#'   similarity scores. One needs to input a similarity score matrix to the \code{sim}
#'   argument for this method. The procedure is the same as \code{"rds"} but before 
#'   the Sinkhorn algorithm is applied, the entries of the random matrix are scaled by 
#'   \code{sim}.
#' \item \code{"convex"} returns the doubly stochastic matrix from the last iteration of running the Frank-
#' Wolfe algorithm with convex relaxation initialized at the barycenter. For this method, one needs to
#' input two graphs \code{A} and \code{B}, as well as \code{seeds} if applicable.
#' }
#' 
#' @examples
#' ss <- matrix(c(5, 4, 4, 3), nrow = 2)
#' # initialize start matrix without soft seeds
#' init_start(start = "bari", nns = 5)
#' init_start(start = "rds", nns = 3)
#' init_start(start = "rds_perm_bari", nns = 5)
#' init_start(start = "rds_from_sim", nns = 3, sim = matrix(runif(9), 3))
#'
#' # initialize start matrix with soft seeds
#' init_start(start = "bari", nns = 5, ns = 1, soft_seeds = ss)
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
#' \donttest{
#' # FW graph matching with incorrect seeds to start at convex start
#' init_start(start = "convex", nns = 8, ns = 2, soft_seeds = ss, A = g1, B = g2, seeds = seeds)
#' }
#'
#' @export
init_start <- function(start, nns, ns = 0, soft_seeds = NULL, seeds = NULL, ...){

  soft_seeds <- check_seeds(soft_seeds, nns + ns)$seeds
  nss <- nrow(soft_seeds)
  
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
        stop(e, "\nNote: functions passed to init_start must have",
          " at least the arguments nns, ns, and softs_seeds")
      })

    # if we get back a full size matrix then just return
    if (all(dim(start) == nns)) {
      return(start)
    }
    # otherwise add in soft seeds below after checking size
    if (all(dim(start) != nns - nss)) {
      stop("Functions passed to init start must return",
        " a square matrix-like object with dimension ", nns,
        " or", nns - nss)
    }
  } else if (start == "bari"){
    start <- bari_start(nns, ns, soft_seeds)
  } else if (start == "rds") {
    start <- rds_sinkhorn_start(nns, ns, soft_seeds, ...)
  } else if (start == "rds_perm_bari") {
    start <- rds_perm_bari_start(nns, ns, soft_seeds, ...)
  } else if (start == "rds_from_sim"){
    start <- rds_from_sim_start(nns, ns, soft_seeds, ...)
  } else if (start == "convex") {
    # start at bari with soft seeds
    start <- init_start("bari", nns, ns, soft_seeds, seeds)
    # match and pull out doubly stochastic
    start <- gm(..., seeds = seeds, method = "convex", start = start)$soft

    # don't add back in soft seeds b/c we've used them for convex
    # maybe message/warning about this being a silly thing
    # to do
    if (!is.null(seeds)) {
      # needed to avoid check problems
      nonseeds <- check_seeds(seeds, nv = nns + ns)$nonseeds
      start <- start[nonseeds$A, nonseeds$B]
    }
    return(start)
  } else {
    stop('start must be either a matrix, function, or one of "bari", "rds",
      "rds_perm_bari", "rds_from_sim", "convex"')
  }
  if(ns > 0 & is.null(seeds)) {
    seeds <- seq(ns)
  }
  add_soft_seeds(start, nns, ns, soft_seeds, seeds)
}


# Func that takes a nns - nss dim matrix and returns
# a nns dim matrix that incorporates soft seeds
add_soft_seeds <- function(start, nns, ns, soft_seeds, hard_seeds) {
  if (is.null(soft_seeds))
    return(start)

  hard_seeds <- check_seeds(hard_seeds, nv = nns + ns)$seeds

  reindex <- function(s, hs) s - sum(hs < s)

  cs <- check_seeds(soft_seeds, nv = nns + ns)
  
  seeds_g1 <- as.integer(sapply(cs$seeds$A, reindex, hs = hard_seeds$A))
  seeds_g2 <- as.integer(sapply(cs$seeds$B, reindex, hs = hard_seeds$B))
  
  cs <- check_seeds(cbind(A = seeds_g1, B = seeds_g2), nv = nns)
  nonseeds_g1 <- cs$nonseeds$A
  nonseeds_g2 <- cs$nonseeds$B
  nss <- nrow(cs$seeds)

  reorderA <- order(c(nonseeds_g1, seeds_g1))
  reorderB <- order(c(nonseeds_g2, seeds_g2))

  new_start <- pad(start, nss)[reorderA, reorderB]

  # Hack to avoid message about inefficiently treating single elements
  suppressMessages(
    new_start[cbind(seeds_g1, seeds_g2)] <- 1
  )
  new_start
}

