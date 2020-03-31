
graph_match_sPATH <- function(A, B, seeds = NULL, similarity = NULL, 
                              alpha = 0.5, epsilon = 1, lap_method = NULL){
  A <- A[]
  B <- B[]
  Sym <- Matrix::isSymmetric(A) & Matrix::isSymmetric(B)
  # same cardinality
  totv1 <- dim(A)[1]
  totv2 <- dim(B)[1]
  nv <- max(totv1, totv2)
  
  seed_check <- check_seeds(seeds, nv)
  seeds <- seed_check$seeds
  nonseeds <- seed_check$nonseeds
  nns <- nv - nrow(seeds)
  
  degree_A <- Matrix::colSums(A)
  degree_B <- Matrix::colSums(B)
  D_A <- Matrix::Diagonal(length(degree_A), degree_A)
  D_B <- Matrix::Diagonal(length(degree_B), degree_B)
  L_A <- D_A - A
  L_B <- D_B - B
  delta_cal <- function(x, y){
    (y - x) ^ 2
  }
  delta <- outer(Matrix::diag(D_A), Matrix::diag(D_B), delta_cal)
  
  Asn <- A[seeds$A, nonseeds$A]
  Ans <- A[nonseeds$A, seeds$A]
  Ann <- A[nonseeds$A, nonseeds$A]
  Bsn <- B[seeds$B, nonseeds$B]
  Bns <- B[nonseeds$B, seeds$B]
  Bnn <- B[nonseeds$B, nonseeds$B]
  deltann <- delta[nonseeds$A, nonseeds$B]
  L_Asn <- L_A[seeds$A, nonseeds$A]
  L_Ans <- L_A[nonseeds$A, seeds$A]
  L_Ann <- L_A[nonseeds$A, nonseeds$A]
  L_Bsn <- L_B[seeds$B, nonseeds$B]
  L_Bns <- L_B[nonseeds$B, seeds$B]
  L_Bnn <- L_B[nonseeds$B, nonseeds$B]
  
  lap_method <- set_lap_method(lap_method, totv1, totv2)
  
  # lambda=0: convex objective function
  if(is.null(similarity)){
    convex_m <- graph_match_convex(A, B, seeds = seeds, 
                                   start = "bari", tol = 1e-10)
  } else{
    start_sim <- sinkhorn(similarity[nonseeds$A, nonseeds$B])
    convex_m <- graph_match_convex(A, B, seeds = seeds, #start = start_sim,
                                   similarity = similarity, tol = 1e-10)
  }
  P <- convex_m$P[nonseeds$A, nonseeds$B]
  lambda <- 0
  dlambda <- dlambda_min <-  2e-2 ####### IF ADAPTIVE METHOD APPLIED, NO NEED TO WORRY ABOUT INITIALIZATION
  iter <- 0
  
  while (lambda < 1){
    iter <- iter + 1
    # calculate F_cv & F_cc with current P
    if(Sym){
      AB <- Ans %*% Bsn
      APB <- Ann %*% P %*% Bnn
      F_cv <- 4 * sum(P * AB) + 2 * sum(P * APB)
      LALB <- L_Ans %*% L_Bsn
      LAPLB <- L_Ann %*% P %*% L_Bnn
      F_cc <- sum(t(deltann) * P) + 4 * sum(P * LALB) + 2 * sum(P * LAPLB)
    } else{
      ABsn <- Matrix::t(Asn) %*% Bsn
      ABns <- Ans %*% Matrix::t(Bns)
      APBt <- Ann %*% P %*% Matrix::t(Bnn)
      F_cv <- 2*(sum(P * ABsn) + sum(P * ABns) + sum(P * APBt))
      LALBns <- L_Ans %*% Matrix::t(L_Bns)
      LALBsn <- Matrix::t(L_Asn) %*% L_Bsn
      LAtPLB <- Matrix::t(L_Ann) %*% P %*% L_Bnn
      F_cc <- sum(t(deltann) * P) + 2 * (sum(P * LALBns) + 
                                           sum(P * LALBsn) + sum(P * LAtPLB))
    }
    
    # calculate F_old(P_old)
    if(!is.null(similarity)){
      F_sim <- sum(similarity[nonseeds$A, nonseeds$B] * P)
      F_P <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    } else{
      F_P <- (1 - lambda) * F_cv + lambda * F_cc
    } 
    
    # calculate F_new(P_old)
    lambda <- lambda + dlambda
    if(lambda > 1){
      lambda <- 1
    }
    if(!is.null(similarity)){
      F_P_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    } else{
      F_P_new <- (1 - lambda) * F_cv + lambda * F_cc
    }
    
    # dlambda-adaptive
    while (sum(abs(F_P - F_P_new)) < epsilon && lambda < 1) {
      dlambda <- 2 * dlambda
      lambda <- lambda + dlambda
      if(lambda > 1){
        lambda <- 1
        break
      }
      if(!is.null(similarity)){
        F_P_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
      } else{
        F_P_new <- (1 - lambda) * F_cv + lambda * F_cc
      }    
    }
    
    # Frank-Wolfe
    if(Sym){
      Grad_cv <- 4 * (AB + APB)
      Grad_cc <- Matrix::t(deltann) + 4 * (LALB + LAPLB)
    } else{
      Grad_cv <- 2 * (ABsn + ABns + APBt + Matrix::t(Ann) %*% P %*% Bnn)
      Grad_cc <- Matrix::t(deltann) + 2 * (LALBns + LALBsn + LAtPLB + L_Ann %*% P %*% Matrix::t(L_Bnn))
    }
    Grad <- (1 - lambda) * Grad_cv / max(abs(Grad_cv)) + lambda * Grad_cc / max(abs(Grad_cc))
    if(!is.null(similarity)){
      Grad <- alpha * Grad + (1 - alpha) * similarity[nonseeds$A, nonseeds$B]
    }
    Grad <- Grad - min(Grad)
    ind <- do_lap(Grad, lap_method)
    ind2 <- cbind(1:nns, ind)
    Pdir <- Matrix::Diagonal(nns)[ind, ]
    
    # step size
    if(sum(P != Pdir) != 0){
      delta_P <- P - Pdir
      c1 <- sum(Bsn * (Asn %*% delta_P)) + sum(delta_P * (Ans %*% Matrix::t(Bns))) +
        sum(Pdir * (Ann %*% delta_P %*% Matrix::t(Bnn))) + sum(delta_P * (Ann %*% Pdir %*% Matrix::t(Bnn)))
      d1 <- sum(delta_P * (Ann %*% delta_P %*% Matrix::t(Bnn)))  
      if(Sym){
        c2 <- -sum(Matrix::t(deltann) * delta_P) - 2 * (2 * sum(delta_P * (L_Ans %*% Matrix::t(L_Bns))) + 
                                                          sum(L_Ann * (delta_P %*% L_Bnn %*% Matrix::t(Pdir))) + 
                                                          sum(L_Ann * (Pdir %*% L_Bnn %*% Matrix::t(delta_P))))
      } else{
        c2 <- -sum(Matrix::t(deltann) * delta_P) - 
          2 * (sum(delta_P * (L_Ans %*% Matrix::t(L_Bns))) + sum(delta_P * (Matrix::t(L_Asn) %*% L_Bsn)) +
                 sum(L_Ann * (delta_P %*% L_Bnn %*% Matrix::t(Pdir))) + sum(L_Ann * (Pdir %*% L_Bnn %*% Matrix::t(delta_P))))
      }
      d2 <- sum(L_Ann * (delta_P %*% L_Bnn %*% Matrix::t(delta_P)))
      scale <- max(c(c1, c2, d1, d2))
      c1 <- c1/scale
      d1 <- d1/scale
      c2 <- c2/scale
      d2 <- d2/scale
      
      if(is.null(similarity)){
        a <- lambda * c2 - 2 * (1 - lambda) * c1
        b <- 4 * lambda * d2 + 4 * (1 - lambda) * d1
      } else{
        a <- lambda * alpha * c2 - 2 * c1 * alpha * (1-lambda) + 
          (1 - alpha) * sum(similarity[nonseeds$A, nonseeds$B] * delta_P)
        b <- 4 * d1 * alpha * (1 - lambda) + 4 * d2 * lambda * alpha
      }
      if(a == 0 && b == 0){
        beta <- 0
      } else{
        beta <- a/b
      }
      if(beta > 1){
        beta <- 1
      } else if(beta < 0){
        beta <- 0
      }  
      P <- beta * P + (1 - beta) * Pdir
    } 
  }
  
  corr_ns <- do_lap(P, lap_method)
  
  D_ns <- P
  corr_ns <- do_lap(P, lap_method)
  corr <- 1:nv
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B
  P <- Matrix::Diagonal(nv)[corr, ]
  D <- P
  D[nonseeds$A, nonseeds$B] <- D_ns
  
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = nrow(seeds), 
            P = P, D = D, iter = iter)
  z
}
