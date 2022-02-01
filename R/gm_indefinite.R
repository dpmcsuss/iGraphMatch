#' @title Frank-Wolfe Graph Matching Methods
#'
#' @description Match two given graphs, returns a list of graph matching
#'   results, including matching correspondence vector of \eqn{G_2} with respect
#'   to \eqn{G_1}, doubly stochastic matrix and permutation matrix.
#'
#' @param A A matrix, igraph object, or list of either.
#' @param B A matrix, igraph object, or list of either.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix or a data frame, with the first
#'   column being the indices of \eqn{G_1} and the second column being the
#'   corresponding indices of \eqn{G_2}.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#'   character value like "bari", "rds" or "convex" to initialize the starting matrix.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similarities.
#' @param tol A number. Tolerance of edge disagreements.
#' @param max_iter A number. Maximum number of replacing matches.
#' @param lap_method Choice for lap method. One of "lapjv", "lapmod", or "clue".
#'
#' @rdname gm_fw
#'
#' @return \code{graph_match_indefinite}, \code{graph_match_convex},
#'   \code{graph_match_sinkhorn} and \code{graph_match_PATH}.
#'   return an object of class "\code{\link{graphMatch}}" which is a list containing the following
#'   components:
#'
#'   \describe{
#'     \item{corr_A}{matching correspondence in \eqn{G_1}}
#'     \item{corr_B}{matching correspondence in \eqn{G_2}}
#'     \item{soft}{the doubly stochastic matrix from the last iteration with which one can
#'           extract more than one matching candidates}
#'     \item{iter}{number of iterations until convergence or reaches the \code{max_iter}}
#'     \item{max_iter}{Maximum number of replacing matches}
#'     \item{lap_method}{Choice for solving the LAP}
#'     \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
#'   }
#'
#' @references V. Lyzinski and D. E. Fishkind and M. Fiori and J. T. Vogelstein and C. E. Priebe
#' and G. Sapiro (2016), \emph{Graph Matching: Relax at Your Own Risk}. IEEE TPAMI, pages 60-73.
#' @references V. Lyzinski and D. E. Fishkind and C. E. Priebe (2014), \emph{Seeded Graph Matching
#' for Correlated Erdos-Renyi Graphs}.J. Mach. Learn. Res., pages 3513-3540.
#'
#'
#' @examples
#'
#' # match G_1 & G_2 with some known node pairs as seeds
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:10 <= 3
#' GM_bari <- gm(g1, g2, seeds, method = "indefinite", start = "bari")
#' GM_bari
#' GM_bari[!GM_bari$seeds] # matching correspondence for non-seeds
#'
#' summary(GM_bari, g1, g2, true_label = 1:10)
#'
#' # match G_1 & G_2 with some incorrect seeds
#' hard_seeds <- matrix(c(4,6,5,4),2)
#' seeds <- rbind(as.matrix(check_seeds(seeds, nv = 10)$seeds),hard_seeds)
#' GM_badseed <- gm(g1, g2, seeds, method = "indefinite")
#'
#' GM_badseed[] # get the corresponding permutation matrix
#' GM_badseed %*% g2 # permute the second graph according to match result: PBP^T
#' GM_badseed$soft # doubly stochastic matrix from the last step of Frank-Wolfe iterations
#' GM_badseed$iter # number of iterations
#' GM_badseed$max_iter # preset maximum number of iterations: 20
#'
#' # match two multi-layer graphs
#' gp_list <- replicate(3, sample_correlated_gnp_pair(20, .3, .5), simplify = FALSE)
#' A <- lapply(gp_list, function(gp)gp[[1]])
#' B <- lapply(gp_list, function(gp)gp[[2]])
#'
#' match_multi_layer <- gm(A, B, seeds = 1:10, method = "indefinite", start = "bari", max_iter = 20)
#' summary(match_multi_layer, A, B)
#'
#' @keywords internal
graph_match_indefinite <- function(A, B, seeds = NULL,
  similarity = NULL, start = "bari",
  max_iter = 20, lap_method = NULL) {

  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  nonseeds <- check_seeds(seeds, c(totv1, totv2))$nonseeds
  ns <- nrow(seeds)
  nn <- totv1 - ns
  
  # TODO: remove [[1]]
  # TODO: modify initstart to take nn1 and nn2
  P <- init_start(start = start, nns = totv2 - ns, ns = ns,
    A = A[[1]], B = B[[1]], seeds = seeds)

  P <- P[seq(totv1 - ns), seq(totv2 - ns)]

  iter <- 0
  toggle <- TRUE

  # make a random permutation
  rp <- sample(totv2 - ns)
  rpmat <- Matrix::Diagonal(totv2 - ns)[rp, ]

  # similarity and seed to non-seed
  similarity <- similarity[seq(totv1 - ns), seq(totv2 - ns)] %*% Matrix::t(rpmat) +
    get_s_to_ns(A, B, seeds, nonseeds, rp)

  # keep only nonseeds
  A <- A[nonseeds$A, nonseeds$A]
  B <- B[nonseeds$B, nonseeds$B][rp, rp]

  lap_method <- set_lap_method(lap_method, totv1, totv2)

  while(toggle && iter < max_iter){
    iter <- iter + 1
    # non-seed to non-seed info
    tAnn_P_Bnn <- ml_sum(crossprod(A, P %*% B))

    Grad <- tAnn_P_Bnn + similarity +
      ml_sum(tcrossprod(A %*% P, B))

    ind <- do_lap(Grad, lap_method)[seq(totv1 - ns)]

    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(totv2 - ns)
    Pdir <- Pdir[ind, ]
    tAnn_Pdir_Bnn <- ml_sum(crossprod(A, B[ind, ]))

    c <- innerproduct(tAnn_P_Bnn, P)
    d <- innerproduct(tAnn_Pdir_Bnn, P) + sum(tAnn_P_Bnn[ind2])
    e <- sum(tAnn_Pdir_Bnn[ind2])
    u <- innerproduct(P, similarity)
    v <- sum(similarity[ind2])
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


  if(iter == max_iter){
    warning("Frank-Wolfe iterations reached the maximum iteration, convergence may not occur.")
  }

  corr_ns <- do_lap(P, lap_method)[seq(totv1 - ns)]


  # undo rand perm here
  corr_ns <- rp[corr_ns]

  corr <- 1:totv1
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B

  reorderA <- order(c(nonseeds$A, seeds$A))
  reorderB <- order(c(nonseeds$B, seeds$B))

  D <- pad(P %*% rpmat, ns)[reorderA, reorderB]
  if (is(D, "splrMatrix")) {
    D@x[seeds$A, seeds$B] <- Matrix::Diagonal(ns)
  } else {
    D[seeds$A, seeds$B] <- Matrix::Diagonal(ns)
  }
  cl <- match.call()

  graphMatch(
    corr = data.frame(corr_A = 1:totv1, corr_B = corr),
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
