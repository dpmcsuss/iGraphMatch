#' @rdname graph_match_methods
#' @return \code{graph_match_soft_percolation} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B}, the number of seeds and the order of nodes getting matched.
#'
#' @examples
#' # match G_1 & G_2 using soft percolation graph matching method
#' seeds <- 1:5
#' graph_match_soft_percolation(g1, g2, seeds, r = 2, max_iter = 2)
#'
#' @export
#'
#'
graph_match_soft_percolation <- function(A, B, seeds, r = 2, max_iter = 100){

  # this will make the graphs be matrices if they are igraph objects
  if(igraph::is.igraph(A)){
    weighted <- igraph::is.weighted(A)
  } else{
    if(min(A) < 0){
      weighted <- TRUE
    } else{
      weighted <- max(A) > 1
    }
  }
  A <- A[]
  B <- B[]
  
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  n <- max(totv1, totv2)
  P <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
  seeds <- check_seeds(seeds, n)$seeds
  ns <- nrow(seeds)
  seeds_ori <- seeds
  P[as.matrix(seeds)] <- 1
  
  # initialization of score matrix M & MM
  if(weighted){
    M <- Matrix::Matrix(0, totv1, totv2)
    for(i in 1:nrow(seeds)){
      A_adj <- which(A[seeds$A[i],]>0)
      B_adj <- which(B[seeds$B[i],]>0)
      if(length(A_adj) != 0 && length(B_adj) != 0){
        mark <- outer(A[seeds$A[i],A_adj], B[seeds$B[i],B_adj], cal_mark)
        M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
      }
    }
  } else{
    M <- (Matrix::t(A) %*% P %*% B + A %*% P %*% Matrix::t(B)) / 2
  }
  M[seeds$A,] <- -n
  M[,seeds$B] <- -n
  MM <- M # score matrix w. socres to matched pairs be -infinity

  # initialization for checking replacing cycle
  num <- 0 # # of removal seeds
  remove <- rbind(c(0,0),c(0,0)) # list of removed seeds
  remove_by <- rbind(c(0,0),c(0,0))
  colnames(remove) <- paste0(c("A","B"))
  cyc <- FALSE

  # percolate
  Z <- seeds # set of matched pairs
  ZZ <- c(0,0) # set of dominant conflict existing matches
  while(max(MM)>=r & num<=max_iter){
    # locate best match
    max_ind <- which(MM==max(MM), arr.ind = TRUE)
    conflict_log <- conflict_check(Z, max_ind, logical = TRUE)
    sum_conf <- sum(conflict_log)
    if(sum_conf>0 && sum_conf<length(conflict_log)){
      max_ind <- max_ind[!conflict_log,] # subset of non-conflict matches
    }
    if(!is.null(nrow(max_ind))){
      rnum <- sample(nrow(max_ind),1)
      max_ind <- max_ind[rnum,] # solve tie: give priority to non-conflict matches
    }

    # non-conflict new match
    if(sum_conf != length(conflict_log)){
      Z <- rbind(Z,max_ind)

      # correct MM caused by ZZ
      if(!is.null(nrow(ZZ))){
        MM[ZZ$A,] <- M[ZZ$A,]
        MM[,ZZ$B] <- M[,ZZ$B]
        ZZ <- c(0,0)
      }

      # update mark matrix M & MM
      if(weighted){
        A_adj <- which(A[max_ind[1],]>0)
        B_adj <- which(B[max_ind[2],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
          MM[A_adj, B_adj] <- MM[A_adj, B_adj] + mark
        }
      } else{
        Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
        Pi[max_ind[1], max_ind[2]] <- 1 
        delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
        M <- M + delta
        MM <- MM + delta
      }
      MM[max_ind[1], max_ind[2]] <- -n

    } else{ # conflict new match: only when all the ties correspond to conflict match

      conf_row_ind <- conflict_check(Z, matrix(max_ind,1), logical = FALSE)
      conf_ind <- Z[conf_row_ind,]
      if(nrow(conf_ind==2)==2){ # conflict with two existing matches
        score1 <- M[conf_ind[1,1], conf_ind[1,2]]
        score2 <- M[conf_ind[2,1], conf_ind[2,2]]
        score <- max(score1,score2)
      } else{
        score <- M[conf_ind$A,conf_ind$B]
      }

      if(M[max_ind[1], max_ind[2]]>score){ #replace

        num <- num + 1

        # check cycle
        if(length(conf_row_ind)==2){
          remove <- rbind(remove,Z[conf_row_ind[1],])
          cyc_remove <- check_cycle(remove, Z[conf_row_ind[2],])
          remove <- rbind(remove, Z[conf_row_ind[2],])
          cyc_remove_by <- check_cycle(remove_by,max_ind)
          remove_by <- rbind(remove_by, max_ind)
        } else{
          cyc_remove <- check_cycle(remove, Z[conf_row_ind,])
          remove <- rbind(remove, Z[conf_row_ind,])
          cyc_remove_by <- check_cycle(remove_by,max_ind)
          remove_by <- rbind(remove_by, max_ind)
        }
        cyc <- cyc_remove & cyc_remove_by

        Z <- Z[-conf_row_ind,] # remove conflict match
        Z <- rbind(Z,max_ind) # add new match

        # correct MM caused by ZZ
        if(!is.null(nrow(ZZ))){
          MM[ZZ$A,] <- M[ZZ$A,]
          MM[,ZZ$B] <- M[,ZZ$B]
          ZZ <- c(0,0)
        }

        # update mark matrix: subtract removed seed's effect
        if(weighted){
          for (i in 1:length(conf_row_ind)) {
            A_adj <- which(A[Z$A[conf_row_ind[i]],]>0)
            B_adj <- which(B[Z$B[conf_row_ind[i]],]>0)
            if(length(A_adj) != 0 && length(B_adj) != 0){
              mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
              M[A_adj, B_adj] <- M[A_adj, B_adj] - mark
              MM[A_adj, B_adj] <- MM[A_adj, B_adj] - mark
              MM[conf_ind[i,1], conf_ind[i,2]] <- M[conf_ind[i,1], conf_ind[i,2]]
            }
          }
        } else{
          Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
          Pi[as.matrix(conf_ind)] <- 1 
          delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
          M <- M - delta
          MM <- MM - delta
          MM[as.matrix(conf_ind)] <- M[as.matrix(conf_ind)]
        }

        # update mark matrix M & MM: add new match's effect
        if(weighted){
          A_adj <- which(A[max_ind[1],]>0)
          B_adj <- which(B[max_ind[2],]>0)
          if(length(A_adj) != 0 && length(B_adj) != 0){
            mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
            M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
            MM[A_adj, B_adj] <- MM[A_adj, B_adj] + mark
          }
        } else{
          Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
          Pi[max_ind[1], max_ind[2]] <- 1 
          delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
          M <- M + delta
          MM <- MM + delta
        }
        MM[max_ind[1], max_ind[2]] <- -n

      } else{ # choose another qualified match
        ZZ <- rbind(ZZ, conf_ind)
        MM[conf_ind[,1],] <- -n
        MM[,conf_ind[,2]] <- -n
      }
    }

  }# end while: percolate

  # matching result
  order <- order(Z$A)
  corr <- Z[order,]
  names(corr) <- c("corr_A","corr_B")
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = ns, order = order)
  z
}

conflict_check <- function(Matches, ind, logical = TRUE){

  if(logical == TRUE){
    rconflict <- ind[,1] %in% Matches$A
    cconflict <- ind[,2] %in% Matches$B
    conflict <- rconflict | cconflict
  } else{
    rconflict <- ind[1] == Matches$A
    rind <- which(rconflict==TRUE)
    cconflict <- ind[2] == Matches$B
    cind <- which(cconflict==TRUE)
    conflict <- c(rind,cind)
  }
  conflict
}

check_cycle <- function(rem, new){
  row <- which(rem[,1]==unlist(new[1]))
  col <- which(rem[,2]==unlist(new[2]))
  occ <- intersect(row,col)
  occ <- c(occ,length(rem[,1])+1)
  if(length(occ)==3){
    if(occ[3]-occ[2]==occ[2]-occ[1]){
      cycle1 <- rem[occ[1]:(occ[2]-1),]
      cycle2 <- rem[occ[2]:(occ[3]-1),]
      if(sum(cycle1[,1]==cycle2[,1])+sum(cycle1[,2]==cycle2[,2])==occ[3]-occ[1]){
        result <- TRUE
      } else{
        result <- FALSE
      }
    } else{
      result <- FALSE
    }
  } else{
    result <- FALSE
  }

  result
}