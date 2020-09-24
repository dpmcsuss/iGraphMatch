cal_mark <- function(x,y){
  1 - abs(x - y) / max(x, y)
}


#' @title Percolation Graph Matching Methods
#' @rdname gm_perco
#' 
#' @return \code{graph_match_percolation} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B}, the number of seeds and the order of nodes getting matched.
#' 
#'
#' @param A A matrix, igraph object, or list of either.
#' @param B A matrix, igraph object, or list of either. 
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param r A number. Threshold of neighboring pair scores.
#'
#' @references L. Yartseva and M. Grossglauser (2013), \emph{On the performance
#'   of percolation graph matching}. COSN, Boston, MA, USA, pages 119–130.
#'
#' @examples
#' # match G_1 & G_2 using percolation graph matching method
#' seeds <- 1:5
#' graph_match_percolation(g1, g2, seeds, r = 2)
#'
#' @export
#'
#'
graph_match_percolation <- function (A, B, seeds, r = 2)
{
  graph_pair <- check_graph(A, B, same_order = FALSE, as_list = FALSE)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  
  directed <- !(isSymmetric(A) && isSymmetric(B))

  n <- max(totv1, totv2)
  seeds <- check_seeds(seeds, n)$seeds
  ns <- nrow(seeds)
  Z <- seeds #matched nodes
  M <- matrix(0, totv1, totv2) #marks matrix
  
  # mark neighbors
  for(i in 1:nrow(seeds)){
    A_adj <- which(A[seeds$A[i],]>0)
    B_adj <- which(B[seeds$B[i],]>0)
    if(length(A_adj) != 0 && length(B_adj) != 0){
      mark <- outer(A[seeds$A[i],A_adj], B[seeds$B[i], B_adj], cal_mark)
      M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
    }
    if(directed){
      A <- Matrix::t(A)
      B <- Matrix::t(B)
      A_adj <- which(A[seeds$A[i],]>0)
      B_adj <- which(B[seeds$B[i],]>0)
      if(length(A_adj) != 0 && length(B_adj) != 0){
        mark <- outer(A[seeds$A[i],A_adj], B[seeds$B[i], B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      }
      A <- Matrix::t(A)
      B <- Matrix::t(B)
    }
  }
  M[seeds$A,] <- -n
  M[,seeds$B] <- -n
  
  while(max(M) >= r){
    max_ind <- Matrix::which(M == max(M), arr.ind = TRUE)
    max_ind <- max_ind[sample(nrow(max_ind),1),]
    
    Z <- rbind(Z,max_ind)
    
    # update mark matrix
    A_adj <- which(A[max_ind[1],]>0)
    B_adj <- which(B[max_ind[2],]>0)
    mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
    M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
    if(directed){
      A <- Matrix::t(A)
      B <- Matrix::t(B)
      A_adj <- which(A[max_ind[1],]>0)
      B_adj <- which(B[max_ind[2],]>0)
      mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
      M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      A <- Matrix::t(A)
      B <- Matrix::t(B)
    }
    
    M[max_ind[1],] <- -n
    M[,max_ind[2]] <- -n
  }
  
  order <- order(Z$A)
  corr <- Z[order,]
  names(corr) <- c("corr_A","corr_B")
  rownames(corr) <- paste0(as.character(1:nrow(corr)))
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = nrow(seeds), order = order)
  z
}