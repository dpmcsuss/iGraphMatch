cal_mark <- function(x,y){
  1 - abs(x - y) / max(abs(x), abs(y))
}


#' @title Percolation Graph Matching Methods
#' @rdname gm_perco
#'
#' @return \code{graph_match_percolation} returns an object of class "gm" which is a
#'   list containing the following components:
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
#' @param ExpandWhenStuck A logical. TRUE if expand the seed set when Percolation algorithm
#'   stops before matching all the vertices.
#'
#' @references L. Yartseva and M. Grossglauser (2013), \emph{On the performance
#'   of percolation graph matching}. COSN, Boston, MA, USA, pages 119–130.
#' @references E. Kazemi, S. H. Hassani, and M. Grossglauser (2015),
#' \emph{Growing a graph matching from a handful of seeds}. Proc. of the VLDB
#' Endowment, 8(10):1010–1021.
#'
#' @examples
#' # match G_1 & G_2 using percolation graph matching method
#' cgnp_pair <- sample_correlated_gnp_pair(n = 20, corr =  0.5, p =  0.8)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:10 <= 3
#' GM_perco <- gm(g1, g2, seeds, method = "percolation", r = 2, ExpandWhenStuck = FALSE)
#' GM_perco
#'
#' # matching accuracy with the true alignment being the identity
#' mean(GM_perco$corr_A == GM_perco$corr_B)
#' GM_perco$match_order
#'
#' summary(GM_perco, g1, g2, true_label = 1:20)
#' plot(g1[], g2[], GM_perco)
#'
#' # expand when stuck
#' GM_exp <- gm(g1, g2, seeds, method = "percolation", r = 4, ExpandWhenStuck = TRUE)
#' GM_exp
#'
#'
#'
graph_match_percolation <- function (A, B, seeds,
                            similarity = NULL, r = 2,
                            ExpandWhenStuck = FALSE) {

  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  n <- max(totv1, totv2)
  ns <- nrow(seeds)
  nn <- n - ns
  nc <- length(A)

  nonseeds <- check_seeds(seeds, n)$nonseeds
  seeds_ori <- seeds
  Z <- seeds #matched nodes

  M <- Matrix(0, n, n) #marks matrix
  if(sum(abs(similarity)) != 0){
    M[nonseeds$A, nonseeds$B] <- similarity
  }
  if(nrow(seeds_ori) != 0){
    M[seeds_ori$A,] <- -n * n
    M[,seeds_ori$B] <- -n * n
  }

  iter <- 1
  while(ns != 0 | (sum(abs(similarity)) != 0 & iter == 1)){

    # mark neighbors
    if(ns != 0){
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
    }


    # choose pairs with marks ge r
    while(max(M) >= r){
      max_ind <- Matrix::which(M == max(M), arr.ind = TRUE)
      if(nrow(max_ind) != 1){
        degree_diff <- 0
        for (ch in 1:nc) {
          degree_diff <- degree_diff + abs(rowSums(A[[ch]])[max_ind[,1]]-rowSums(B[[ch]])[max_ind[,2]])
        }
        max_ind <- max_ind[which(degree_diff == min(degree_diff)),]
        if(is.vector(max_ind) == FALSE){
          max_ind <- max_ind[sample(nrow(max_ind),1),]
        }
      }

      # update mark matrix
      for( ch in 1:nc ){
        directed <- !(isSymmetric(A[[ch]]) && isSymmetric(B[[ch]]))

        A_adj <- which(A[[ch]][max_ind[1],]!=0)
        B_adj <- which(B[[ch]][max_ind[2],]!=0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[[ch]][max_ind[1],A_adj], B[[ch]][max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
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

      M[max_ind[1],] <- -n * n
      M[,max_ind[2]] <- -n * n
      max_ind <- data.frame(A = max_ind[1], B = max_ind[2])
      Z <- rbind(Z, max_ind)
    }

    # deferred percolation gm
    iter <- iter + 1
    ns <- 0
    if(ExpandWhenStuck){
      seeds_old <- seeds
      seeds <- which(M > 0 & M < r, arr.ind = TRUE)
      seeds <- data.frame(A = seeds[,1], B = seeds[,2])
      ns <- nrow(seeds)

      if( (nrow(seeds) == nrow(seeds_old)) && (sum(seeds == seeds_old) == 2*nrow(seeds))){
        break
      }
    }

  }


  # matching result
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
      seeds = seeds_ori
    )
  )
}
