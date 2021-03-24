#' @title Frank-Wolfe Graph Matching Methods
#'
#' @description Match two given graphs, returns a list of graph matching
#'   results, including matching correspondence vector of \eqn{G_2} with respect
#'   to \eqn{G_1}, doubly stochastic matrix and permutation matrix.
#'
#' @param A A matrix, 'igraph' object, or list of either.
#' @param B A matrix, 'igraph' object, or list of either.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix or a data frame, with the first
#'   column being the indices of \eqn{G_1} and the second column being the
#'   corresponding indices of \eqn{G_2}.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#'   character value like "bari" or "convex" to initialize the starting matrix.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similarities.
#' @param tol A number. Tolerance of edge disagreements.
#' @param max_iter A number. Maximum number of replacing matches equals to
#'   max_iter times number of total vertices of \eqn{G_1}.
#' @param lap_method Choice for lap method.
#'
#' @rdname gm_fw
#'
#' @return \code{graph_match_FW}, \code{graph_match_convex} and \code{graph_match
#'   _PATH} return a list of graph matching results, including the graph matching
#'   formula, a data frame containing the matching correspondence between \eqn{G_1}
#'   and \eqn{G_2} named \code{corr_A} and \code{corr_B}, the doubly stochastic
#'   matrix from the last iteration and the permutation matrix after projection,
#'   seeds and number of iterations.
#'
#' @examples
#'
#' # match G_1 & G_2 with some known node pairs as seeds
#' seeds <- 1:10 <= 3
#' graph_match_FW(g1, g2, seeds, start = "bari")
#'
#' # match G_1 & G_2 with some incorrect seeds
#' hard_seeds <- matrix(c(4,6,5,4),2)
#' seeds <- rbind(as.matrix(check_seeds(seeds, nv = 10)$seeds),hard_seeds)
#' graph_match_FW(g1, g2, seeds, start = "convex")
#'
#'  gp_list <- replicate(3, sample_correlated_gnp_pair(20, .3, .5), simplify = FALSE)
#'  A <- lapply(gp_list, function(gp)gp[[1]])
#'  B <- lapply(gp_list, function(gp)gp[[2]])
#'  match <- graph_match_FW(A, B, seeds = 1:10, start = "bari", max_iter = 20)
#'  match$corr
#'
#' @export
#'
graph_match_FW <- function(A, B, seeds = NULL,
  similarity = NULL, start = "bari",
  max_iter = 20, lap_method = NULL) {


  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2

  nv <- nrow(A[[1]])

  seed_check <- check_seeds(seeds, nv)
  seeds <- seed_check$seeds
  nonseeds <- seed_check$nonseeds

  ns <- nrow(seeds)
  nn <- nv - ns

  P <- init_start(start = start, nns = nn, ns = ns,
    A = A[[1]], B = B[[1]], seeds = seeds)

  iter <- 0
  toggle <- TRUE

  # make a random permutation
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]

  # seed to non-seed info
  s_to_ns <- get_s_to_ns(A, B, seeds, nonseeds, rp)
  P <- P[, rp]

  zero_mat <- Matrix::Matrix(0, nn, nn)

  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  similarity <- similarity %*% Matrix::t(rpmat)

  # keep only nonseeds
  A <- A[nonseeds$A, nonseeds$A]
  B <- B[nonseeds$B, nonseeds$B][rp, rp]
  nc <- length(A)

  lap_method <- set_lap_method(lap_method, totv1, totv2)



  while(toggle && iter < max_iter){
    iter <- iter + 1
    # non-seed to non-seed info
    tAnn_P_Bnn <- zero_mat
    for( ch in 1:nc ){
      tAnn_P_Bnn <- tAnn_P_Bnn +
        Matrix::t(A[[ch]]) %*% P %*% B[[ch]]
    }

    Grad <- s_to_ns + tAnn_P_Bnn + similarity
    for(ch in 1:nc){
      Grad <- Grad + A[[ch]] %*% P %*% Matrix::t(B[[ch]])
    }

    ind <- do_lap(Grad, lap_method)

    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(nn)
    Pdir <- Pdir[ind, ]
    ns_Pdir_ns <- zero_mat
    for(ch in 1:nc){
      ns_Pdir_ns <- ns_Pdir_ns +
        Matrix::t(A[[ch]])[, order(ind)] %*% B[[ch]]

    }
    c <- innerproduct(tAnn_P_Bnn, P)
    d <- innerproduct(ns_Pdir_ns, P) + sum(tAnn_P_Bnn[ind2])
    e <- sum(ns_Pdir_ns[ind2])
    u <- innerproduct(P, s_to_ns + similarity)
    v <- sum((s_to_ns + similarity)[ind2])
    if (c - d + e == 0 && d - 2 * e + u - v == 0) {
      alpha <- 0
    } else {
      alpha <- -(d - 2 * e + u - v)/(2 * (c - d + e))
    }
    f0 <- 0
    f1 <- c - e + u - v
    falpha <- (c - d + e) * alpha^2 + (d - 2 * e + u - v) *
      alpha

    if (alpha < 1 && alpha > 0 &&
        falpha > f0 && falpha > f1) {
      P <- alpha * P + (1 - alpha) * Pdir
    } else if (f0 > f1) {
      P <- Pdir
    } else {
      toggle <- F
    }
  }

  corr_ns <- do_lap(P, lap_method)


  # undo rand perm here
  corr_ns <- rp[corr_ns]

  corr <- 1:nv
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B

  reorderA <- order(c(nonseeds$A, seeds$A))
  reorderB <- order(c(nonseeds$B, seeds$B))

  D <- pad(P %*% rpmat, ns)[reorderA, reorderB]
  if (is(D, "splrMatrix")) {
    D@x[seeds$A, seeds$B] <- Matrix::Diagonal(ns)
     # <- P[seeds$A, seeds$B]
  } else {
    D[seeds$A, seeds$B] <- Matrix::Diagonal(ns)
     # <- P[seeds$A, seeds$B]
  }
  cl <- match.call()

  graphMatch(
    corr = data.frame(corr_A = 1:nv, corr_B = corr),
    nnodes = c(totv1, totv2),
    call = cl,
    detail = list(
      iter = iter,
      max_iter = max_iter,
      lap_method = lap_method,
      seeds = seeds,
      soft = D
    )
  )
}

#' @rdname gm_fw
#' @export
gm_indefinite <- graph_match_FW

