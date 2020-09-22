#' @rdname graph_match_methods
#' 
#' @return \code{graph_match_PATH} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B}, the number of seeds if any, the permutation matrix and the
#'   doubly stochastic matrix before projection onto the permutation set. 
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
graph_match_PATH <- function(A, B, similarity = NULL, seeds = NULL, alpha = .5, epsilon = 1){
  graph_pair <- check_graph(A, B, as_list = FALSE)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  
  D_A <- Matrix::Diagonal(length(rowSums(A)), rowSums(A))
  D_B <- Matrix::Diagonal(length(rowSums(B)), rowSums(B))
  A <- A[]
  B <- B[]
  L_A <- D_A - A
  L_B <- D_B - B
  n <- nrow(A)

  seed_check <- check_seeds(seeds, n)
  seeds <- seed_check$seeds
  nonseeds <- seed_check$nonseeds


  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  
  # alpha=0, convex relaxation
  convex_m <- graph_match_convex(A, B, similarity = similarity, seeds = seeds, tol = 1e-10)
  P <- convex_m$P
  lambda <- 0
  dlambda <- dlambda_min <-  1e-5
  #toggle <- TRUE
  delta_cal <- function(x, y){
    (y - x) ^ 2
  }
  delta <- outer(diag(D_A), diag(D_B), delta_cal)
  iter <- 0
  
  lap_method <- set_lap_method(NULL, totv1, totv2)

  while (lambda < 1) {
    iter <- iter + 1
    # dlambda-adaptation
    F_cv <- (Matrix::norm(A %*% P - P %*% B, type = "F")) ^ 2
    L <- Matrix::kronecker(Matrix::t(L_B), Matrix::t(L_A))
    F_cc <- - sum(t(delta) %*% P) - 2 * t(Matrix::c.sparseVector(P)) %*% 
      L %*% Matrix::c.sparseVector(P)
  
    F_sim <- sum(similarity * P)
    F <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    
    lambda <- lambda + dlambda
    if(!is.null(similarity)){
      F_sim <- sum(similarity * P)
      F_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    } else{
      F_new <- (1 - lambda) * F_cv + lambda * F_cc
    }
    
    while (sum(abs(F - F_new)) < epsilon && lambda < 1) {
      dlambda <- 2 * dlambda
      lambda <- lambda + dlambda
      if(lambda > 1){
        lambda <- 1
        break
      }
      F_sim <- sum(similarity * P)
      F_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    }
    while (sum(abs(F - F_new)) > epsilon && dlambda != dlambda_min) {
      if(lambda > 1){
        lambda <- 1
        break
      } else{
        dlambda <- dlambda / 2
        lambda <- lambda - dlambda
      }
      F_sim <- sum(similarity * P)
      F_new <- alpha * ((1 - lambda) * F_cv + lambda * F_cc) + (1 - alpha) * F_sim
    }
    
    # Frank-Wolfe 
    AtA <- t(A) %*% A
    BBt <- B %*% t(B)
    tA_P_B <- Matrix::t(A) %*% P %*% B
    Grad_cv <- 2 * (AtA %*% P + P %*% BBt - tA_P_B - A %*% P %*% t(B))
    Grad_cc <- - t(delta) - 2 * Matrix::t(L_A) %*% P %*% L_B
    Grad <- (1 - lambda) * Grad_cv + lambda * Grad_cc
    Grad <- alpha * Grad + (1 - alpha) * similarity

    ind <- do_lap(Grad, lap_method)
    ind2 <- cbind(1:n, ind)
    Pdir <- Matrix::Diagonal(n)[ind, ]
    
    delta_P <- P - Pdir
    C <- A %*% delta_P - delta_P %*% B
    D <- A %*% Pdir - Pdir %*% B
    aq <- innerproduct(C, C)
    bq <- innerproduct(C, D)
    vec_delta_P <- Matrix::c.sparseVector(delta_P)
    vec_Pdir <- Matrix::c.sparseVector(Pdir)
    c <- sum(t(delta) * delta_P)
    # NOT SURE WHAT SHOULD GO HERE
    # WAS F BEFORE
    e <- Matrix::t(vec_Pdir) %*% L %*% vec_Pdir
    u <- Matrix::t(vec_Pdir) %*% L %*% vec_delta_P
    v <- Matrix::t(vec_delta_P) %*% L %*% vec_delta_P
    a <- 2 * (lambda - 1) * bq + lambda * (c - e + u)
    b <- 2 * (1 - lambda) * aq - 4 * lambda * v
    if(a[1,1] == 0 && b[1,1] == 0){
      alpha <- 0
    } else{
      alpha <- (2 * (lambda - 1) * bq + lambda * (c - e + u)) / 
        (2 * (1 - lambda) * aq - 4 * lambda * v)
      alpha <- alpha[1,1]
    }
    if(alpha > 1){
      alpha <- 1
    } else if(alpha < 0){
      alpha <- 0
    }
    P <- alpha * P + (1 - alpha) * Pdir
  }
  
  D <- P
  corr <- do_lap(P, lap_method)
  P <- Matrix::Diagonal(n)[corr,]
  
  if(!is.null(seeds)){
    ns <- nrow(check_seeds(seeds, n)$seeds)
  } else{
    ns <- 0
  }
  cl <- match.call()
  z <- list(call = cl, corr = data.frame(corr_A = 1:nrow(A), corr_B = corr), ns = ns, 
            P = P, D = D, iter = iter, lambda = lambda)
  z
}