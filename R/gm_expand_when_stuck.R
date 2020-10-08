cal_mark <- function(x,y){
  1 - abs(x - y) / max(x, y)
}


#' @rdname gm_perco
#' @return \code{graph_match_ExpandWhenStuck} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B}, the number of seeds and the order of nodes getting matched.
#'
#' @references E. Kazemi, S. H. Hassani, and M. Grossglauser (2015),
#' \emph{Growing a graph matching from a handful of seeds}. Proc. of the VLDB
#' Endowment, 8(10):1010–1021.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 using Expand When Stuck graph matching method
#' seeds <- 1:5
#' graph_match_ExpandWhenStuck(g1, g2, seeds, r = 2)
#'
#' @export
#'
#'
graph_match_ExpandWhenStuck <- function(A, B, seeds, 
                                        similarity = NULL, r = 2){

  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  nc <- length(A)

  n <- max(totv1, totv2)
  seeds <- check_seeds(seeds, nv = max(totv1, totv2))
  nonseeds <- seeds$nonseeds
  seeds_ori <- seeds <- seeds$seeds
  ns <- nrow(seeds)
  Z <- seeds #matched nodes
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  similarity <- (similarity - min(similarity)) / (max(similarity) - min(similarity))
  if(is.na(similarity)[1,1]){
    M <- matrix(0, totv1, totv2) #marks matrix
  } else{
    M <- similarity
  }
  M[seeds_ori$A,] <- -n * n
  M[,seeds_ori$B] <- -n * n

  # deferred percolation graph matching
  while(nrow(seeds) != 0){
    # mark neighbors
    for(ch in 1:nc){
      directed <- !(isSymmetric(A[[ch]]) && isSymmetric(B[[ch]]))
      
      for(i in 1:nrow(seeds)){
        A_adj <- which(A[[ch]][seeds$A[i],]>0)
        B_adj <- which(B[[ch]][seeds$B[i],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[[ch]][seeds$A[i],A_adj], B[[ch]][seeds$B[i], B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
        if(directed){
          A[[ch]] <- Matrix::t(A[[ch]])
          B[[ch]] <- Matrix::t(B[[ch]])
          A_adj <- which(A[[ch]][seeds$A[i],]>0)
          B_adj <- which(B[[ch]][seeds$B[i],]>0)
          if(length(A_adj) != 0 && length(B_adj) != 0){
            mark <- outer(A[[ch]][seeds$A[i],A_adj], B[[ch]][seeds$B[i], B_adj], cal_mark)
            M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
          }
          A[[ch]] <- Matrix::t(A[[ch]])
          B[[ch]] <- Matrix::t(B[[ch]])
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
        A_adj <- which(A[[ch]][max_ind[1],]>0)
        B_adj <- which(B[[ch]][max_ind[2],]>0)
        mark <- outer(A[[ch]][max_ind[1],A_adj], B[[ch]][max_ind[2],B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        if(directed){
          A[[ch]] <- Matrix::t(A[[ch]])
          B[[ch]] <- Matrix::t(B[[ch]])
          A_adj <- which(A[[ch]][max_ind[1],]>0)
          B_adj <- which(B[[ch]][max_ind[2],]>0)
          mark <- outer(A[[ch]][max_ind[1],A_adj], B[[ch]][max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
          A[[ch]] <- Matrix::t(A[[ch]])
          B[[ch]] <- Matrix::t(B[[ch]])
        }
      }
  
      M[max_ind[1],] <- -n * n
      M[,max_ind[2]] <- -n * n
      max_ind <- data.frame(A = max_ind[1], B = max_ind[2])
      Z <- rbind(Z, max_ind)
    }

    seeds <- which(M > 0 & M < r, arr.ind = TRUE)
    seeds <- data.frame(A = seeds[,1], B = seeds[,2])

  }

  # matching result
  order <- order(Z$A)
  corr <- Z[order(Z$A),]
  names(corr) <- c("corr_A","corr_B")
  rownames(corr) <- paste0(as.character(1:nrow(corr)))
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = ns, order = order)
  z
}
