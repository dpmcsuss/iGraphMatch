#' @title Spectral Graph Matching Methods: Umeyama Algorithm
#' @rdname gm_Umeyama
#'
#' @param A A matrix, 'igraph' object, or list of either.
#' @param B A matrix, 'igraph' object, or list of either.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similarities.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#'
#' @return \code{graph_match_Umeyama} returns an object of class "gm" which is a list
#'   containing the following components:
#'
#'   \describe{
#'     \item{corr_A}{matching correspondence in \eqn{G_1}}
#'     \item{corr_B}{matching correspondence in \eqn{G_2}}
#'     \item{soft}{the functional similarity score matrix with which one can extract
#'       more than one matching candidates}
#'     \item{lap_method}{Choice for solving the LAP}
#'     \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
#'   }
#'
#' @references S. Umeyama (1988), \emph{An eigendecomposition approach to weighted
#'   graph matching problems}. IEEE TPAMI. USA, pages 695-703.
#'
#' @examples
#' # match G_1 & G_2 using Umeyama algorithm
#' G <- sample_correlated_gnp_pair(10, .9, .5)
#' g1 <- G$graph1
#' g2 <- G$graph2
#' startm <- matrix(0, 10, 10)
#' diag(startm)[1:4] <- 1
#'
#' GM_Umeyama <- gm(g1, g2, similarity = startm, method = "Umeyama")
#' GM_Umeyama
#' # generate the corresponding permutation matrix
#' GM_Umeyama[]
#'
#' summary(GM_Umeyama, g1, g2)
#' # visualize the edge-wise matching performance
#' plot(g1, g2, GM_Umeyama)
#' plot(g1[], g2[], GM_Umeyama)
#'
#'
#' @keywords internal
graph_match_Umeyama <- function(A, B, seeds = NULL,
                                similarity = NULL){

  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  nv <- max(totv1, totv2)
  nonseeds <- check_seeds(seeds, nv)$nonseeds
  ns <- nrow(seeds)
  nn <- nv - ns
  nc <- length(A)
  Grad <- 0

  for( ch in 1:nc ){
    if(!isSymmetric(as.matrix(A[[ch]])) | !isSymmetric(as.matrix(B[[ch]]))){
      # construct Hermitian matrices by adjacency matrices
      A[[ch]] <- as.matrix((A[[ch]] + Matrix::t(A[[ch]]))/2) + as.matrix((A[[ch]] - Matrix::t(A[[ch]]))/2)*1i
      B[[ch]] <- as.matrix((B[[ch]] + Matrix::t(B[[ch]]))/2) + as.matrix((B[[ch]] - Matrix::t(B[[ch]]))/2)*1i
    }

    U_A <- eigen(A[[ch]])$vectors
    U_B <- eigen(B[[ch]])$vectors
    AB <- Matrix::tcrossprod(abs(U_B), abs(U_A))
    Grad <- Grad + AB[nonseeds$A, nonseeds$B]
  }


  Grad <- Grad + Matrix::t(similarity)
  Grad <- Grad - min(Grad)

  # make a random permutation
  nn <- nrow(A[[1]]) - nrow(seeds)
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]
  Grad <- Grad[,rp]

  lap_method <- set_lap_method(NULL, totv1, totv2)
  ind <- do_lap(Grad, lap_method)

  # undo rand perm here
  ind <- rp[ind]
  corr <- data.frame(corr_A = c(seeds$A, nonseeds$A), corr_B = c(seeds$B, nonseeds$B[ind]))
  corr <- corr[order(corr$corr_A),]
  cl <- match.call()

  graphMatch(
    corr = corr,
    nnodes = c(totv1, totv2),
    call = cl,
    detail = list(
      lap_method = lap_method,
      seeds = seeds,
      soft = Grad
    )
  )
}
