

get_graph_triple <- function(g, weight, first_graph){
  if(first_graph){
    w <- sign(weight) * w
    list(w[1] * g[], w[2] * g[],
      splr_sparse_plus_constant(- w[3] * g[], w[3]))
  }
  else{
    list(g[], splr_sparse_plus_constant(- g[], 1), g[])
  }
}


graph_match_FW_multi_reward <- function(A, B, weight, ...){
  if( is.igraph(A) ){
    A <- list(A)
  }
  if( is.igraph(B) ){
    B <- list(B)
  }
  if(!is.list(A) && !is.list(B)){
    A <- list(A)
    B <- list(B)
  }
  A <- unlist(lapply(A, get_graph_triple,
    weight = weight, first_graph = TRUE), recursive = FALSE)
  B <- unlist(lapply(B, get_graph_triple,
    weight = weight, first_graph = FALSE), recursive = FALSE)

  graph_match_FW_multi(A, B, ...)
}

              

graph_match_percolation_multi <- function (A, B, start = NULL, similarity = NULL, r = 2, alpha = 10) 
{
  A <- lapply(A, function(Al) Al[])
  B <- lapply(B, function(Bl) Bl[])
  
  n_A <- ncol(A[[1]])
  n_B <- ncol(B[[1]])
  P <- start * alpha
  Z <- c(0,0)
  nc <- length(A)
  M <- similarity
  for(ch in 1:nc){
    M <- M + t(A[[ch]]) %*% P %*% B[[ch]] + A[[ch]] %*% P %*% t(B[[ch]])
    gc()
  }
  
  while (ifelse(is.na(max(M)), 0, max(M)) >= r + max(similarity)) {
    max_ind <- which(M == max(M), arr.ind = TRUE)
    max_ind <- max_ind[sample(nrow(max_ind), 1), ]
    
    P[max_ind[1],] <- 0
    P[,max_ind[2]] <- 0
    P[max_ind[1], max_ind[2]] <- 1
    P[P>0] <- 1
    P <- diag(1 / rowSums(P)) %*% P
    
    M <- similarity
    for(ch in 1:nc){
      M <- M + t(A[[ch]]) %*% P %*% B[[ch]] + A[[ch]] %*% P %*% t(B[[ch]])
      gc()
    }
    #if(min(B)<0){
    #  delta_norm <- 2*delta_norm - matrix(1, nrow = n_A, ncol = n_B)
    #}
    Z <- rbind(Z, max_ind)
    M[Z[,1], ] <- -n_B
    M[,Z[,2]] <- -n_B
  }
  
  corr <- Z[-1,]
  if(is.null(nrow(corr))){
    corr <- Matrix(corr, ncol = 2)
  }
  colnames(corr) <- c("corr_A", "corr_B")
  corr
}





graph_match_ExpandWhenStuck_multi <- function(A, B, start = NULL, similarity = NULL, r = 2, alpha = 5){
  # this will make the graphs be matrices if they are igraph objects
  A <- lapply(A, function(Al) Al[])
  B <- lapply(B, function(Bl) Bl[])
  
  n_A <- ncol(A[[1]])
  n_B <- ncol(B[[1]])
  P <- start * alpha
  Z <- c(0,0)
  nc <- length(A)
  M <- similarity
  for(ch in 1:nc){
    M <- M + t(A[[ch]]) %*% P %*% B[[ch]] + A[[ch]] %*% P %*% t(B[[ch]])
    gc()
  }
  
  seeds <- Matrix(0,2,2)
  # deferred percolation graph matching
  while(nrow(seeds) != 0){
    if(sum(seeds[1,])!=0){
      Pi <- Matrix(0, nrow = n_A, ncol = n_B)
      Pi[seeds] <- 1
      for(ch in 1:nc){
        M <- M + t(A[[ch]]) %*% Pi %*% B[[ch]] + A[[ch]] %*% Pi %*% t(B[[ch]])
        gc()
      }
    }
    
    # choose pairs with marks ge r
    while(max(M) >= r + max(similarity)){
      max_ind <- which(M == max(M), arr.ind = TRUE)
      max_ind <- max_ind[sample(nrow(max_ind), 1), ]
      
      # update mark matrix
      P[max_ind[1],] <- 0
      P[,max_ind[2]] <- 0
      P[max_ind[1], max_ind[2]] <- 1
      P[P>0] <- 1
      P <- diag(1 / rowSums(P)) %*% P
      if(sum(seeds[1,])!=0){
        P[seeds] <- 1 
      }
      
      M <- similarity
      for(ch in 1:nc){
        M <- M + t(A[[ch]]) %*% P %*% B[[ch]] + A[[ch]] %*% P %*% t(B[[ch]])
        gc()
      }
      #if(min(B)<0){
      #  delta_norm <- 2*delta_norm - matrix(1, nrow = n_A, ncol = n_B)
      #}
      Z <- rbind(Z, max_ind)
      M[Z[,1], ] <- -n_B
      M[,Z[,2]] <- -n_B
    }
    
    seeds_old <- seeds
    seeds <- which(M > 1e5+1, arr.ind = TRUE)
    
    if(nrow(seeds) == nrow(seeds_old)){
      if(sum(seeds == seeds_old)==2*nrow(seeds)){
        break
      }
    }
    if(nrow(Z)==n_A + 1){
      break
    }
  }
  
  # matching result
  corr <- Z[-1,]
  if(is.null(nrow(corr))){
    corr <- Matrix(corr, ncol = 2)
  }
  colnames(corr) <- c("corr_A","corr_B")
  corr
}


graph_match_mutual_multi <- function(A, B, start = NULL, similarity = NULL, alpha = 0.2, max_iter = 50){
  # this will make the graphs be matrices if they are igraph objects
  A <- lapply(A, function(Al) Al[])
  B <- lapply(B, function(Bl) Bl[])
  
  match <- Matrix(0,2,2)
  mutual3 <- match
  mutual21 <- match
  mutual22 <- match
  a0 <- 0 # number of seeds
  iter <- 0
  
  n_A <- ncol(A[[1]])
  n_B <- ncol(B[[1]])
  nc <- length(A)
  n <- min(n_A, n_B)
  
  if(alpha==0){
    while(iter < max_iter & nrow(match) < n){
      P <- start
      P[match[,1],] <- 0
      P[,match[,2]] <- 0
      P[match] <- 1 # update permutation matrix by current matches
      P <- diag(1 / rowSums(P)) %*% P
      M <- similarity
      for(ch in 1:nc){
        M <- M + t(A[[ch]]) %*% P %*% B[[ch]] + A[[ch]] %*% P %*% t(B[[ch]])
        gc()
      }
      col_ind <- max.col(M[(a0+1):n_A,(a0+1):n_B], ties.method = "random") + a0 # strategy for solving ties: random
      row_ind <- max.col(t(M[(a0+1):n_A,(a0+1):n_B]), ties.method = "random") + a0 # strategy for solving ties: random
      col_max <- cbind(1:n_A, col_ind)
      row_max <- cbind(row_ind, 1:n_B)
      candidate <- rbind(col_max, row_max)
      dup <- duplicated(candidate)
      match <- candidate[dup,]
      iter <- iter + 1
      
    }
  } else{
    while(iter < max_iter & nrow(match) < n){
      P <- start
      P[match[,1],] <- 0
      P[,match[,2]] <- 0
      P[as.matrix(match)] <- 1 # update permutation matrix by current matches
      P <- diag(1 / rowSums(P)) %*% P
      P[as.matrix(mutual3)] <- alpha
      P[as.matrix(mutual21)] <- alpha * 3 / 2
      P[as.matrix(mutual22)] <- alpha * 3 / 2
      M <- similarity
      for(ch in 1:nc){
        M <- M + t(A[[ch]]) %*% P %*% B[[ch]] + A[[ch]] %*% P %*% t(B[[ch]])
        gc()
      }
      
      col_ind <- max.col(as.matrix(M[(a0+1):n_A,(a0+1):n_B]), ties.method = "random") + a0 # strategy for solving ties: random
      row_ind <- max.col(as.matrix(t(M[(a0+1):n_A,(a0+1):n_B])), ties.method = "random") + a0 # strategy for solving ties: random
      row_max1 <- cbind(1:n_A, col_ind) # pairs of nodes with highest score in each row
      col_max1 <- cbind(row_ind, 1:n_B) # pairs of nodes with highest score in each column
      
      MM <- M
      MM[row_max1] <- 0
      col_ind2 <- max.col(MM[(a0+1):n_A,(a0+1):n_B], ties.method = "random") + a0
      row_max2 <- cbind(1:n_A, col_ind2)
      
      MM <- M
      MM[col_max1] <- 0
      row_ind2 <- max.col(t(MM[(a0+1):n_A,(a0+1):n_B]), ties.method = "random") + a0
      col_max2 <- cbind(row_ind2, 1:n_B) 
      
      candidate1 <- rbind(col_max1, row_max1)
      dup1 <- duplicated(candidate1)
      match <- candidate1[dup1,]
      
      candidate21 <- rbind(col_max1, row_max2)
      dup21 <- duplicated(candidate21)
      mutual21 <- candidate21[dup21,]
      
      candidate22 <- rbind(col_max2, row_max1)
      dup22 <- duplicated(candidate22)
      mutual22 <- candidate22[dup22,]
      
      candidate3 <- rbind(col_max2, row_max2)
      dup3 <- duplicated(candidate3)
      mutual3 <- candidate3[dup3,]
      
      iter <- iter + 1
      
    }
  }
  
  corr <- match
  if(is.null(nrow(corr))){
    corr <- Matrix(corr, ncol = 2)
  }
  colnames(corr) <- c("corr_A", "corr_B") 
  corr
}
