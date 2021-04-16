cal_mark <- function(x,y){
  1 - abs(x - y) / max(abs(x), abs(y))
}


#' @title Percolation Graph Matching Methods
#' @rdname gm_perco
#'
#' @return \code{graph_match_percolation} and \code{graph_match_ExpandWhenStuck}
#'   return an object of class "gm" which is a list containing the following
#'   components:
#'
#'   \describe{
#'     \item{corr_A}{matching correspondence in \eqn{G_1}}
#'     \item{corr_B}{matching correspondence in \eqn{G_2}}
#'     \item{match_order}{the order of vertices getting matched}
#'     \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
#'   }
#'
#'
#'
#'
#' @param A A matrix, 'igraph' object, or list of either.
#' @param B A matrix, 'igraph' object, or list of either.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similarities.
#' @param r A number. Threshold of neighboring pair scores.
#'
#' @references L. Yartseva and M. Grossglauser (2013), \emph{On the performance
#'   of percolation graph matching}. COSN, Boston, MA, USA, pages 119–130.
#'
#' @examples
#' # match G_1 & G_2 using percolation graph matching method
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:10 <= 3
#' GM_perco <- graph_match_percolation(g1, g2, seeds, r = 2)
#' GM_perco
#'
#' # matching accuracy with the true alignment being the identity
#' mean(GM_perco$corr_A == GM_perco$corr_B)
#' GM_perco$match_order
#'
#' summary(GM_perco, g1, g2)
#' plot(g1[], g2[], GM_perco)
#'
#' @export
#'
#'
graph_match_percolation <- function (A, B, seeds,
                                     similarity = NULL, r = 2) {

  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  nc <- length(A)

  n <- max(totv1, totv2)
  seeds <- check_seeds(seeds, nv = max(totv1, totv2))
  nonseeds <- seeds$nonseeds
  seeds <- seeds$seeds
  ns <- nrow(seeds)
  if(ns == 0 & is.null(similarity)){
    stop("at least one of seeds and similarity score should be known for this method.")
  }
  Z <- seeds #matched nodes
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  M <- Matrix(0, n, n) #marks matrix
  if(sum(abs(similarity)) != 0){
    M[nonseeds$A, nonseeds$B] <- similarity
  }


  # mark neighbors
  if(nrow(seeds) != 0){
    for(ch in 1:nc){
      directed <- !(isSymmetric(A[[ch]]) && isSymmetric(B[[ch]]))

      for(i in 1:nrow(seeds)){
        A_adj <- which(A[[ch]][seeds$A[i],]!=0)
        B_adj <- which(B[[ch]][seeds$B[i],]!=0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[[ch]][seeds$A[i],A_adj], B[[ch]][seeds$B[i], B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
        if(directed){
          A[[ch]] <- Matrix::t(A[[ch]])
          B[[ch]] <- Matrix::t(B[[ch]])
          A_adj <- which(A[[ch]][seeds$A[i],]!=0)
          B_adj <- which(B[[ch]][seeds$B[i],]!=0)
          if(length(A_adj) != 0 && length(B_adj) != 0){
            mark <- outer(A[[ch]][seeds$A[i],A_adj], B[[ch]][seeds$B[i], B_adj], cal_mark)
            M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
          }
          A[[ch]] <- Matrix::t(A[[ch]])
          B[[ch]] <- Matrix::t(B[[ch]])
        }
      }
    }

    M[seeds$A,] <- -n
    M[,seeds$B] <- -n
  }


  while(max(M) >= r){
    max_ind <- Matrix::which(M == max(M), arr.ind = TRUE)
    max_ind <- max_ind[sample(nrow(max_ind),1),]

    Z <- rbind(Z,max_ind)

    # update mark matrix
    for( ch in 1:nc ){
      A_adj <- which(A[[ch]][max_ind[1],]!=0)
      B_adj <- which(B[[ch]][max_ind[2],]!=0)
      if(length(A_adj) != 0 && length(B_adj) != 0){
        mark <- outer(A[[ch]][max_ind[1],A_adj], B[[ch]][max_ind[2],B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      }
      directed <- !(isSymmetric(A[[ch]]) && isSymmetric(B[[ch]]))
      if(directed){
        A[[ch]] <- Matrix::t(A[[ch]])
        B[[ch]] <- Matrix::t(B[[ch]])
        A_adj <- which(A[[ch]][max_ind[1],]!=0)
        B_adj <- which(B[[ch]][max_ind[2],]!=0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[[ch]][max_ind[1],A_adj], B[[ch]][max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
        A[[ch]] <- Matrix::t(A[[ch]])
        B[[ch]] <- Matrix::t(B[[ch]])
      }
    }

    M[max_ind[1],] <- -n
    M[,max_ind[2]] <- -n
  }

  order <- order(Z[,1])
  corr <- Z[order,]
  names(corr) <- c("corr_A","corr_B")
  rownames(corr) <- paste0(as.character(1:nrow(corr)))

  cl <- match.call()
  graphMatch(
    corr = corr,
    nnodes = c(totv1, totv2),
    call = cl,
    detail = list(
      match_order = order,
      seeds = seeds
    )
  )
}
