#' @title Spectral Graph Matching Methods: IsoRank Algorithm
#' @rdname gm_isorank
#'
#' @param A A matrix, igraph object, or list of either.
#' @param B A matrix, igraph object, or list of either.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similarities.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param max_iter A number. Maximum number of replacing matches.
#' @param lap_method Choice of method to extract mapping from score matrix.
#'   One of "greedy" or "LAP".
#'
#' @return \code{graph_match_IsoRank} returns an object of class "\code{\link{graphMatch}}" which is a list
#'   containing the following components:
#'
#'   \describe{
#'     \item{corr_A}{matching correspondence in \eqn{G_1}}
#'     \item{corr_B}{matching correspondence in \eqn{G_2}}
#'     \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
#'     \item{soft}{the functional similarity score matrix obtained from the power method
#'       with which one can extract more than one matching candidates}
#'     \item{match_order}{the order of vertices getting matched}
#'     \item{lap_method}{Method for extracting node mapping}
#'   }
#'
#'
#' @references R. Singh, J. Xu, B. Berger (2008), \emph{Global alignment of
#' multiple protein interaction networks with application to functional
#' orthology detection}. Proceedings of the National Academy of Science. USA, pages 12763-12768.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 using IsoRank algorithm
#' startm <- as.matrix(init_start(start = "bari", nns = 10, soft_seeds = 1:4))
#'
#' GM_IsoRank <- gm(g1, g2, similarity = startm, method = "IsoRank", lap_method = "greedy")
#' GM_IsoRank
#' summary(GM_IsoRank, g1, g2, true_label = 1:10)
#'
#' GM_IsoRank[] # get the corresponding permutation matrix
#' GM_IsoRank %*% g2 # permute the second graph according to match result: PBP^T
#' GM_IsoRank %*% g2[] # output permuted matrix
#'
#' # Visualize the edge-wise matching performance
#' plot(g1, g2, GM_IsoRank)
#' plot(g1[], g2[], GM_IsoRank)
#'
#'
#' @keywords internal
graph_match_IsoRank <- function(A, B, seeds = NULL, similarity,
                                max_iter = 50, lap_method = "greedy"){

  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  nv <- max(totv1, totv2)
  nonseeds <- check_seeds(seeds, nv)$nonseeds
  ns <- nrow(seeds)
  nn <- nv - ns
  nc <- length(A)

  R <- E <- similarity / sum(abs(similarity))
  tol <- 1e-2
  R_tot <- Matrix(0, nrow(R), ncol(R))

  for( ch in 1:nc ) {

    iter <- 1
    diff <- 1

    # computing transition matrix A
    colS_A <- Matrix::colSums(A[[ch]])
    colS_B <- Matrix::colSums(B[[ch]])
    A[[ch]] <- A[[ch]] %*% Matrix::Diagonal(nrow(A[[ch]]), ifelse(colS_A == 0, 0, 1/colS_A))
    B[[ch]] <- B[[ch]] %*% Matrix::Diagonal(nrow(B[[ch]]), ifelse(colS_B == 0, 0, 1/colS_B))

    # computing R by power method
    while(diff > tol & iter <= max_iter){

      AR <- A[[ch]] %*% R %*% Matrix::t(B[[ch]])
      AR <- AR + E
      R_new <- AR / sum(abs(AR))
      diff <- sum(abs(R-R_new))
      iter <- iter + 1
      R <- R_new
    }

    R_tot <- R_tot + R
  }

  # find GNA
  R <- R_tot[nonseeds$A, nonseeds$B]
  if(lap_method == "greedy"){
    corr <- NULL
    while (max(R)>0) {
      max_ind <- Matrix::which(R == max(R), arr.ind = TRUE)
      max_ind <- max_ind[sample(nrow(max_ind), 1), ]
      corr <- rbind(corr, max_ind)
      R[max_ind[1],] <- -1
      R[,max_ind[2]] <- -1
    }
    corr <- data.frame(corr_A = c(seeds$A, nonseeds$A[corr[,1]]),
                       corr_B = c(seeds$B, nonseeds$B[corr[,2]]))
    order <- order(corr$corr_A)
    corr <- corr[order,]
    names(corr) <- c("corr_A","corr_B")
    rownames(corr) <- paste0(as.character(1:nrow(corr)))
    cl <- match.call()

    D <- Matrix(0, nrow(R_tot), ncol(R_tot))
    D[seeds$A, seeds$B] <- diag(nrow(seeds))
    D[nonseeds$A, nonseeds$B] <- R_tot[nonseeds$A, nonseeds$B]

    m <- graphMatch(
      corr = corr,
      nnodes = c(totv1, totv2),
      call = cl,
      detail = list(
        lap_method = lap_method,
        match_order = order,
        seeds = seeds,
        soft = D
      )
    )
  } else if(lap_method == "LAP") {
    # make a random permutation
    nn <- nrow(A[[1]]) - nrow(seeds)
    rp <- sample(nn)
    rpmat <- Matrix::Diagonal(nn)[rp, ]
    R <- R %*% Matrix::t(rpmat)
    # Hungarian alg.
    lap_method <- set_lap_method(NULL, totv1, totv2)
    corr <- do_lap(R - min(R), lap_method)
    # undo rand perm here
    corr <- rp[corr]
    corr <- data.frame(corr_A = c(seeds$A, nonseeds$A), corr_B = c(seeds$B, nonseeds$B[corr]))
    corr <- corr[order(corr$corr_A),]
    names(corr) <- c("corr_A","corr_B")
    rownames(corr) <- paste0(as.character(1:nrow(corr)))
    cl <- match.call()
    D <- Matrix(0, nrow(R_tot), ncol(R_tot))
    D[seeds$A, seeds$B] <- diag(nrow(seeds))
    D[nonseeds$A, nonseeds$B] <- R_tot[nonseeds$A, nonseeds$B]

    m <- graphMatch(
      corr = corr,
      nnodes = c(totv1, totv2),
      call = cl,
      detail = list(
        lap_method = lap_method,
        seeds = seeds,
        soft = D
      )
    )
  }
  return(m)
}
