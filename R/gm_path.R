delta_cal <- function(x, y){
  (y - x) ^ 2
}
#' @rdname gm_fw
#' 
#' @references M. Zaslavskiy, F. Bach and J. Vert (2009), \emph{A Path following
#' algorithm for the graph matching problem}. IEEE Trans Pattern Anal Mach Intell,
#' pages 2227-2242.
#'
#' @param epsilon A small number
#' 
#' @examples
#' # match G_1 & G_2 using PATH algorithm
#' graph_match_PATH(g1, g2)
#'
#' @export
#'
#'
graph_match_PATH <- function(A, B, seeds = NULL, similarity = NULL, 
                             epsilon = 1, tol = 1e-05, max_iter = 20, 
                             lap_method = NULL){
  
  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  nc <- length(A)
  n <- nrow(A)
  
  seeds <- check_seeds(seeds, n)
  nonseeds <- seeds$nonseeds
  seeds <- seeds$seeds
  ns <- nrow(seeds)
  nn <- n - ns
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  
  # lambda=0, convex relaxation
  convex_m <- graph_match_convex(A, B, similarity = similarity, seeds = seeds, 
                                 tol = tol, max_iter = max_iter)
  P <- convex_m[][nonseeds$A, nonseeds$B]
  
  lambda <- 0
  dlambda <- dlambda_min <-  1e-2
  iter <- 0
  lap_method <- set_lap_method(NULL, totv1, totv2)
  
  Grad <- Matrix(0, n, n)
  Sym <- norm <- list()
  Asn <- Ans <- Ann <- list()
  Bsn <- Bns <- Bnn <- list()
  L_Asn <- L_Ans <- L_Ann <- list()
  L_Bsn <- L_Bns <- L_Bnn <- list()
  deltann <- list()
  
  # make a random permutation
  nn <- nrow(A[[1]]) - nrow(seeds)
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]
  similarity <- similarity %*% Matrix::t(rpmat)
  P <- P %*% Matrix::t(rpmat)
  
  ###### multi-layer starts from here
  for (ch in 1:nc) {
    Sym[[ch]] <- isSymmetric(A[[ch]]) && isSymmetric(B[[ch]])
    norm[[ch]] <- sqrt(Matrix::norm(A[[ch]], "F")^2 + Matrix::norm(B[[ch]], "F")^2)
    
    D_A <- Matrix::Diagonal(length(rowSums(A[[ch]])), rowSums(A[[ch]]))
    D_B <- Matrix::Diagonal(length(rowSums(B[[ch]])), rowSums(B[[ch]]))
    L_A <- D_A - A[[ch]]
    L_B <- D_B - B[[ch]]
    delta <- outer(diag(D_A), diag(D_B), delta_cal)
    
    
    Asn[[ch]] <- A[[ch]][seeds$A, nonseeds$A]
    Ans[[ch]] <- A[[ch]][nonseeds$A, seeds$A]
    Ann[[ch]] <- A[[ch]][nonseeds$A, nonseeds$A]
    Bsn[[ch]] <- B[[ch]][seeds$B, nonseeds$B] %*% Matrix::t(rpmat)
    Bns[[ch]] <- rpmat %*% B[[ch]][nonseeds$B, seeds$B]
    Bnn[[ch]] <- rpmat %*% B[[ch]][nonseeds$B, nonseeds$B] %*% Matrix::t(rpmat)
    deltann[[ch]] <- rpmat %*% delta[nonseeds$B, nonseeds$A]
    L_Asn[[ch]] <- L_A[seeds$A, nonseeds$A]
    L_Ans[[ch]] <- L_A[nonseeds$A, seeds$A]
    L_Ann[[ch]] <- L_A[nonseeds$A, nonseeds$A]
    L_Bsn[[ch]] <- L_B[seeds$B, nonseeds$B] %*% Matrix::t(rpmat)
    L_Bns[[ch]] <- rpmat %*% L_B[nonseeds$B, seeds$B]
    L_Bnn[[ch]] <- rpmat %*% L_B[nonseeds$B, nonseeds$B] %*% Matrix::t(rpmat)
  }
  
  
  while (lambda < 1){
    iter <- iter + 1
    
    # calculate F_cv & F_cc with current P
    F_cc <- F_cv <- 0
    Grad_cc <- Grad_cv <- Matrix::Matrix(0, nn, nn)
    
    for (ch in 1:nc) {
      if(Sym[[ch]]){
        AB <- Ans[[ch]] %*% Bsn[[ch]]
        APB <- Ann[[ch]] %*% P %*% Bnn[[ch]]
        F_cv <- F_cv + 4 * sum(P * AB) + 2 * sum(P * APB)
        LALB <- L_Ans[[ch]] %*% L_Bsn[[ch]]
        LAPLB <- L_Ann[[ch]] %*% P %*% L_Bnn[[ch]]
        F_cc <- F_cc + (sum(t(deltann[[ch]]) * P) + 4 * sum(P * LALB) + 
                          2 * sum(P * LAPLB)) / norm[[ch]]
        
        Grad_cv <- Grad_cv + 4 * (AB + APB)
        Grad_cc <- Grad_cc + (Matrix::t(deltann[[ch]]) + 4 * (LALB + LAPLB)) / norm[[ch]]
      } else{
        ABsn <- Matrix::t(Asn[[ch]]) %*% Bsn[[ch]]
        ABns <- Ans[[ch]] %*% Matrix::t(Bns[[ch]])
        APBt <- Ann[[ch]] %*% P %*% Matrix::t(Bnn[[ch]])
        F_cv <- F_cv + 2*(sum(P * ABsn) + sum(P * ABns) + sum(P * APBt))
        LALBns <- L_Ans[[ch]] %*% Matrix::t(L_Bns[[ch]])
        LALBsn <- Matrix::t(L_Asn[[ch]]) %*% L_Bsn[[ch]]
        LAtPLB <- Matrix::t(L_Ann[[ch]]) %*% P %*% L_Bnn[[ch]]
        F_cc <- F_cc + (sum(t(deltann[[ch]]) * P) + 2 * (sum(P * LALBns) + 
                                                           sum(P * LALBsn) + sum(P * LAtPLB))) / norm[[ch]]
        
        Grad_cv <- Grad_cv + 2 * (ABsn + ABns + APBt + Matrix::t(Ann[[ch]]) %*% P %*% Bnn[[ch]])
        Grad_cc <- Grad_cc + (Matrix::t(deltann[[ch]]) + 
                                2 * (LALBns + LALBsn + LAtPLB + L_Ann[[ch]] %*% P %*% Matrix::t(L_Bnn[[ch]]))) / norm[[ch]]
      }
      
    }
    
    # calculate F_old(P_old)
    F_sim <- sum(similarity * P)
    F_P <- (1 - lambda) * F_cv + lambda * F_cc + F_sim 
    
    # calculate F_new(P_old)
    lambda <- lambda + dlambda
    if(lambda > 1){
      lambda <- 1
    }
    F_P_new <- (1 - lambda) * F_cv + lambda * F_cc + F_sim
    
    # dlambda-adaptive
    while (lambda < 1 && sum(abs(F_P - F_P_new)) < epsilon) {
      dlambda <- 2 * dlambda
      lambda <- lambda + dlambda
      if(lambda > 1){
        lambda <- 1
        break
      }
      F_P_new <- (1 - lambda) * F_cv + lambda * F_cc + F_sim  
    }
    # Frank-Wolfe
    Grad <- (1 - lambda) * Grad_cv  + lambda * Grad_cc + similarity
    Grad <- Grad - min(Grad)
    ind <- do_lap(Grad, lap_method)
    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(nn)[ind, ]
    
    # step size
    if(sum(P != Pdir) != 0){
      delta_P <- P - Pdir
      c1 <- c2 <- d1 <- d2 <- 0
      
      for(ch in 1:nc){
        c1 <- c1 + sum(Bsn[[ch]] * (Asn[[ch]] %*% delta_P)) + 
          sum(delta_P * (Ans[[ch]] %*% Matrix::t(Bns[[ch]]))) +
          sum(Pdir * (Ann[[ch]] %*% delta_P %*% Matrix::t(Bnn[[ch]]))) + 
          sum(delta_P * (Ann[[ch]] %*% Pdir %*% Matrix::t(Bnn[[ch]])))
        d1 <- d1 + sum(delta_P * (Ann[[ch]] %*% delta_P %*% Matrix::t(Bnn[[ch]])))  
        if(Sym[[ch]]){
          c2 <- c2 + -sum(Matrix::t(deltann[[ch]]) * delta_P) - 
            2 * (2 * sum(delta_P * (L_Ans[[ch]] %*% Matrix::t(L_Bns[[ch]]))) + 
                   sum(L_Ann[[ch]] * (delta_P %*% L_Bnn[[ch]] %*% Matrix::t(Pdir))) + 
                   sum(L_Ann[[ch]] * (Pdir %*% L_Bnn[[ch]] %*% Matrix::t(delta_P))))
        } else{
          c2 <- c2 + -sum(Matrix::t(deltann[[ch]]) * delta_P) - 
            2 * (sum(delta_P * (L_Ans[[ch]] %*% Matrix::t(L_Bns[[ch]]))) + 
                   sum(delta_P * (Matrix::t(L_Asn[[ch]]) %*% L_Bsn[[ch]])) +
                   sum(L_Ann[[ch]] * (delta_P %*% L_Bnn[[ch]] %*% Matrix::t(Pdir))) + 
                   sum(L_Ann[[ch]] * (Pdir %*% L_Bnn[[ch]] %*% Matrix::t(delta_P))))
        }
        d2 <- d2 + sum(L_Ann[[ch]] * (delta_P %*% L_Bnn[[ch]] %*% Matrix::t(delta_P)))
      }
      
      a <- lambda * c2 - 2 * c1 * (1-lambda) + sum(similarity * delta_P)
      b <- 4 * d1 * (1 - lambda) + 4 * d2 * lambda 
      if(a == 0 && b == 0){
        alpha <- 0
      } else{
        alpha <- a/b
      }
      if(alpha > 1){
        alpha <- 1
      } else if(alpha < 0){
        alpha <- 0
      } 
    } else{
      alpha <- 1
    }
    P <- alpha * P + (1 - alpha) * Pdir
    
  }


  corr_ns <- do_lap(P, lap_method)
  # undo rand perm here
  corr_ns <- rp[corr_ns]
  
  corr <- 1:n
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B
  # D <- P
  # D[nonseeds$A, nonseeds$B] <- D_ns %*% rpmat
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
    call = cl, 
    corr = data.frame(corr_A = seq(n), corr_B = corr),
    nnodes = c(totv1, totv2),
    detail = list(
      seeds = seeds,
      soft = D, 
      iter = iter
    )
  )
}