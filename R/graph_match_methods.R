#' @title Graph Match Methods
#'
#' @description Match two given graphs, returns a list of graph matching
#'   results, including matching correspondence vector of \eqn{G_2} with respect
#'   to \eqn{G_1}, doubly stochastic matrix and permutation matrix.
#'
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   there is no error in seeds input can be a vector of seed indices in
#'   \eqn{G_1}. Or if there exists error in seeds, input in the form of a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#'   character value like "bari" or "convex" to initialize the starting matrix.
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similaities.
#' @param tol A number. Tolerance of edge disagreements.
#' @param r A number. Threshold of neighboring pair scores.
#' @param max_iter A number. Maximum number of replacing matches equals to
#'   max_iter times number of total vertices of \eqn{G_1}.
#' @param alpha A number betwen 0 and 1. Bigger alpha means putting more importance
#'   on the information in network topology over other information such as
#'   similarity scores
#' @param lap_method Choice for lap method.
#' @param method A character. Choice of method to extract mapping from score matrix,
#'   including greedy method and the Hungarian algorithm.
#'
#' @rdname graph_match_methods
#'   
#' @return \code{graph_match_FW} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B} and the number of seeds. 
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 with no seeds
#' graph_match_FW(g1, g2)
#'
#' # match G_1 & G_2 with some known node pairs as seeds
#' seeds <- 1:10 <= 3
#' graph_match_FW(g1, g2, seeds, start = "bari")
#'
#' # match G_1 & G_2 with some incorrect seeds
#' hard_seeds <- matrix(c(4,6,5,4),2)
#' seeds <- rbind(as.matrix(check_seeds(seeds, nv = 10)$seeds),hard_seeds)
#' graph_match_FW(g1, g2, seeds, start = "convex")
#'
#'  gp_list <- replicate(3, sample_correlated_gnp_pair(100, .3, .5), simplify = FALSE)
#'  A <- lapply(gp_list, function(gp)gp[[1]])
#'  B <- lapply(gp_list, function(gp)gp[[2]])
#'  match <- graph_match_FW(A, B, seeds = 1:10, start = "bari", max_iter = 20)
#'  match$corr
#'
#' @export
#'
graph_match_FW <- function(A, B, seeds = NULL,
  start = "convex", max_iter = 20,
  similarity = NULL, lap_method = NULL) {


  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2

  nv <- nrow(A[[1]])

  seed_check <- check_seeds(seeds, nv)
  seeds <- seed_check$seeds
  nonseeds <- seed_check$nonseeds

  ns <- nrow(seeds)
  nn <- nv - ns

  P <- init_start(start = start, nns = nn, ns = ns,
    A = A[[1]], B = B[[1]], seeds = seeds)

  iter <- 0
  toggle <- TRUE

  # make a random permutation
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]

  # seed to non-seed info
  s_to_ns <- get_s_to_ns(A, B, seeds, nonseeds, rp)
  P <- P[, rp]

  zero_mat <- Matrix::Matrix(0, nn, nn)

  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  similarity <- similarity %*% Matrix::t(rpmat)

  # keep only nonseeds
  A <- A[nonseeds$A, nonseeds$A]
  B <- B[nonseeds$B, nonseeds$B][rp, rp]
  nc <- length(A)

  lap_method <- set_lap_method(lap_method, totv1, totv2)



  while(toggle && iter < max_iter){
    iter <- iter + 1
    # non-seed to non-seed info
    tAnn_P_Bnn <- zero_mat
    for( ch in 1:nc ){
      tAnn_P_Bnn <- tAnn_P_Bnn +
        Matrix::t(A[[ch]]) %*% P %*% B[[ch]]
    }

    Grad <- s_to_ns + tAnn_P_Bnn + similarity
    for(ch in 1:nc){
      Grad <- Grad + A[[ch]] %*% P %*% Matrix::t(B[[ch]])
    }

    ind <- do_lap(Grad, lap_method)

    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(nn)
    Pdir <- Pdir[ind, ]
    ns_Pdir_ns <- zero_mat
    for(ch in 1:nc){
      ns_Pdir_ns <- ns_Pdir_ns +
        Matrix::t(A[[ch]])[, order(ind)] %*% B[[ch]]
      
    }
    c <- innerproduct(tAnn_P_Bnn, P)
    d <- innerproduct(ns_Pdir_ns, P) + sum(tAnn_P_Bnn[ind2])
    e <- sum(ns_Pdir_ns[ind2])
    u <- innerproduct(P, s_to_ns + similarity)
    v <- sum((s_to_ns + similarity)[ind2])
    if (c - d + e == 0 && d - 2 * e + u - v == 0) {
      alpha <- 0
    } else {
      alpha <- -(d - 2 * e + u - v)/(2 * (c - d + e))
    }
    f0 <- 0
    f1 <- c - e + u - v
    falpha <- (c - d + e) * alpha^2 + (d - 2 * e + u - v) *
      alpha

    if (alpha < 1 && alpha > 0 &&
        falpha > f0 && falpha > f1) {
      P <- alpha * P + (1 - alpha) * Pdir
    } else if (f0 > f1) {
      P <- Pdir
    } else {
      toggle <- F
    }
  }

  D_ns <- P

  corr_ns <- do_lap(P, lap_method)


  # undo rand perm here
  corr_ns <- rp[corr_ns]

  corr <- 1:nv
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B
  P <- Matrix::Diagonal(nv)[corr, ]
  D <- P
  D[nonseeds$A, nonseeds$B] <- D_ns %*% rpmat

  cl <- match.call()
  list(
    call = cl, 
    corr = data.frame(corr_A = 1:nv, corr_B = corr),
    ns = ns,
    P = P,
    D = D,
    num_iter = iter)
}

#' @rdname graph_match_methods
gm_indefinite <- graph_match_FW

#' @rdname graph_match_methods
#' @return \code{graph_match_convex} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B} and the number of seeds. 
#'
#' @examples
#' seeds <- 1:10 <= 3
#' graph_match_convex(g1, g2, seeds)
#'
#' hard_seeds <- matrix(c(4,6,5,4),2)
#' seeds <- rbind(as.matrix(check_seeds(seeds, 10)$seeds),hard_seeds)
#' \dontrun{
#' # for some reason this fails in check
#' graph_match_convex(g1, g2, seeds)
#' }
#'
#' @export
#'
#'
graph_match_convex <- function(A, B, seeds = NULL, 
  start = "bari", max_iter = 100, similarity = NULL,
  tol = 1e-5, lap_method = NULL) {
  graph_pair <- check_graph(A, B)
  A <- matrix_list(graph_pair[[1]])
  B <- matrix_list(graph_pair[[2]])

  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  nv <- totv1

  seed_check <- check_seeds(seeds, nv)
  seeds <- seed_check$seeds
  nonseeds <- seed_check$nonseeds

  ns <- nrow(seeds)
  nn <- nv - ns


  # make a random permutation
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]

  Asn <- A[seeds$A,nonseeds$A]
  Ann <- A[nonseeds$A,nonseeds$A]
  Ans <- A[nonseeds$A,seeds$A]

  Bsn <- B[seeds$B,nonseeds$B][,rp]
  Bnn <- B[nonseeds$B,nonseeds$B][rp,rp]
  Bns <- B[nonseeds$B,seeds$B][rp,]
  

  zero_mat <- Matrix::Matrix(0, nn, nn)
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  similarity <- similarity %*% Matrix::t(rpmat)

  tol0 <- 1
  if(identical(start, "convex")){
    stop("Cannot start convex with convex. Try \"bari\" or another option.")
  }
  P <- init_start(start = start, nn, ns)[, rp]
  iter <- 0
  toggle <- TRUE



  AtA <- ml_sum(t(Asn) %*% Asn + t(Ann) %*% Ann)
  BBt <- ml_sum(Bns %*% t(Bns) + Bnn %*% t(Bnn))    

  ABns_sn <- ml_sum(Ans %*% t(Bns) + t(Asn) %*% Bsn)


  f <- innerproduct(Ann %*% P - P%*% Bnn,
    Ann %*% P - P%*% Bnn)
    

  lap_method <- set_lap_method(lap_method, totv1, totv2)
  alpha_seq <- NULL
  while(toggle && iter < max_iter){
    f_old <- f
    iter <- iter + 1
    Grad <- ml_sum(
      AtA %*% P + P %*% BBt - ABns_sn +
      -t(Ann) %*% P %*% Bnn - Ann %*% P %*% t(Bnn) +
      similarity)
   

    corr <- do_lap(Grad, lap_method)
    Pdir <- Matrix::Diagonal(nn)[corr,]


    # C <- rbind(Ann,Asn) %*% (P-Pdir) + t((Pdir-P) %*% cbind(Bnn,Bns))
    Cnn <- Ann %*% (P - Pdir) - (P - Pdir) %*% Bnn
    Dnn <- Ann %*% Pdir - Pdir %*% Bnn

    if(ns > 0) {
      Cns <- - (P - Pdir) %*% Bns
      Csn <- Asn %*% (P - Pdir)

      Dns <- Ans - Pdir %*% Bns
      Dsn <- Asn %*% Pdir - Bsn
    } else {
      Dns <- Dsn <- Cns <- Csn <- 0
    }
    aq <- innerproduct(Cnn, Cnn) +
      innerproduct(Cns, Cns) +
      innerproduct(Csn, Csn)
    bq <- innerproduct(Cnn, Dnn) +
      innerproduct(Cns, Dns) +
      innerproduct(Csn, Dsn)
    aopt <- ifelse(aq == 0 && bq == 0, 0,
      ifelse(-bq / aq > 1, 1, -bq/aq))
    alpha_seq <- c(alpha_seq, aopt)
    P_new <- aopt * P + (1 - aopt) * Pdir
    f <- innerproduct(Ann %*% P_new - P_new %*% Bnn,
      Ann %*% P_new - P_new %*% Bnn)

    f_diff <- abs(f - f_old)
    P_diff <- norm(P - P_new, "f")
    P <- P_new

    toggle <- f_diff > tol && f > tol && P_diff > tol
  }



  D_ns <- P

  corr_ns <- do_lap(P, lap_method)
  # undo rand perm here
  corr_ns <- rp[corr_ns]

  corr <- 1:nv
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B
  P <- Matrix::Diagonal(nv)[corr, ]
  D <- P
  D[nonseeds$A, nonseeds$B] <- D_ns %*% rpmat

  cl <- match.call()
  z <- list(
    call = cl,
    corr = data.frame(corr_A = 1:nv, corr_B = corr),
    ns = ns, 
    P = P,
    D = D,
    num_iter = iter,
    alpha_seq = alpha_seq)  
  z
}


#' 
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
#' @rdname graph_match_methods
#' @return \code{graph_match_percolation} returns a list of graph matching results,
#'   including the graph matching formula, a data frame containing the matching 
#'   correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} and 
#'   \code{corr_B}, the number of seeds and the order of nodes getting matched.
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
cal_mark <- function(x,y){
  1 - abs(x - y) / max(x, y)
}
#' @rdname graph_match_methods
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
#' # match G_1 & G_2 using Expand When Stuck graph matching method
#' seeds <- 1:5
#' graph_match_ExpandWhenStuck(g1, g2, seeds, r = 2)
#'
#' @export
#'
#'
graph_match_ExpandWhenStuck <- function(A, B, seeds, r = 2){
  # this will make the graphs be matrices if they are igraph objects
  graph_pair <- check_graph(A, B, same_order = FALSE, as_list = FALSE)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2

  weighted <- any(A != 0 | A != 1 | B != 0 | B != 1)

  n <- max(totv1, totv2)
  P <- Matrix::Matrix(0, nrow = totv1, ncol = totv2)
  seeds <- check_seeds(seeds, n)$seeds
  ns <- nrow(seeds)
  seeds_ori <- seeds
  P[as.matrix(seeds)] <- 1
  M <- Matrix::Matrix(0, totv1, totv2)
  M[seeds_ori$A,] <- -n * n
  M[,seeds_ori$B] <- -n * n
  Z <- seeds

  # deferred percolation graph matching
  while(nrow(seeds) != 0){
    # mark neighbors
    if(weighted){
      for(i in 1:nrow(seeds)){
        A_adj <- which(A[seeds$A[i],]>0)
        B_adj <- which(B[seeds$B[i],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[seeds$A[i],A_adj], B[seeds$B[i],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
      }
    } else{
      Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
      Pi[as.matrix(seeds)] <- 1
      M <- M + (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
    }

    # choose pairs with marks ge r
    while(max(M) >= r){
      max_ind <- Matrix::which(M == max(M), arr.ind = TRUE)
      if(nrow(max_ind) != 1){
        degree_diff <- abs(rowSums(A)[max_ind[,1]]-rowSums(B)[max_ind[,2]])
        max_ind <- max_ind[which(degree_diff == min(degree_diff)),]
        if(is.vector(max_ind) == FALSE){
          max_ind <- max_ind[sample(nrow(max_ind),1),]
        }
      }

      # update mark matrix
      if(weighted){
        A_adj <- which(A[max_ind[1],]>0)
        B_adj <- which(B[max_ind[2],]>0)
        if(length(A_adj) != 0 && length(B_adj) != 0){
          mark <- outer(A[max_ind[1],A_adj], B[max_ind[2],B_adj], cal_mark)
          M[A_adj, B_adj] <- M[A_adj, B_adj] + mark
        }
      } else{
        Pi <- Matrix::Matrix(0, nrow=totv1, ncol = totv2)
        Pi[max_ind[1], max_ind[2]] <- 1 
        delta <- (Matrix::t(A) %*% Pi %*% B + A %*% Pi %*% Matrix::t(B)) / 2
        M <- M + delta
      }
      M[max_ind[1],] <- -n * n
      M[,max_ind[2]] <- -n * n
      max_ind <- data.frame(A = max_ind[1], B = max_ind[2])
      Z <- rbind(Z, max_ind)
    }

    seeds_old <- seeds
    seeds <- which(M > 0 & M < r, arr.ind = TRUE)
    seeds <- data.frame(A=seeds[,1], B=seeds[,2])

    if(nrow(seeds) == nrow(seeds_old)){
      if(sum(seeds == seeds_old)==2*nrow(seeds)){
        break
      }
    }
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
#' @rdname graph_match_methods
#' @return \code{graph_match_IsoRank} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B} and the number of seeds. If choose the greedy method to
#'   extract mapping, the order of nodes getting matched will also be returned.
#'
#' @references R. Singh, J. Xu, B. Berger (2008), \emph{Global alignment of
#' multiple protein interaction networks with application to functional
#' orthology detection}. Proc Natl Acad Sci. USA, pages 12763-12768.
#'
#' @examples
#' # match G_1 & G_2 using IsoRank algorithm
#' startm <- matrix(0, 10, 10)
#' diag(startm)[1:4] <- 1
#' GM_IsoRank <- graph_match_IsoRank(g1, g2, startm, alpha = .3, method = "greedy")
#'
#' @export
#'
graph_match_IsoRank <- function(A, B, similarity, seeds = NULL, 
                                alpha = .5, max_iter = 50, method = "greedy"){
  A <- A[]
  B <- B[]
  
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  
  # padding if two graphs different sizes
  if(totv1 > totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1 < totv2){
    diff <- totv2 - totv1
    A <- pad(A[], diff)
  }
  
  # computing transition matrix A
  colS_A <- Matrix::colSums(A)
  colS_B <- Matrix::colSums(B)
  A <- A %*% Matrix::Diagonal(nrow(A), ifelse(colS_A == 0, 0, 1/colS_A))
  B <- B %*% Matrix::Diagonal(nrow(B), ifelse(colS_B == 0, 0, 1/colS_B))
  E <- similarity / sum(abs(similarity))
  
  # computing R by power method
  R_new <- E
  tol <- 1e-2
  iter <- 1
  diff <- 1
  while(diff > tol & iter <= max_iter){
    
    R <- R_new
    if(alpha>0){
      AR <- A %*% R %*% Matrix::t(B)
      AR <- alpha * AR + (1-alpha) * E
    } else{
      AR <- A %*% R %*% Matrix::t(B)
    }
    R_new <- AR / sum(abs(AR))
    diff <- sum(abs(R-R_new))
    iter <- iter + 1
  }

  seeds_log <- check_seeds(seeds, nv = max(totv1, totv2), logical = TRUE)
  seeds <- check_seeds(seeds, nv = max(totv1, totv2))
  nonseeds <- seeds$nonseeds
  seeds <- seeds$seeds
  R <- R[!seeds_log, !seeds_log]
  R <- as.matrix(R)
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  # find GNA
  if(method == "greedy"){
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
    z <- list(call = cl, corr = corr, ns = nrow(seeds), order = order)
    z
  } else if(method == "LAP"){
    # Hungarian alg.
    lap_method <- set_lap_method(NULL, totv1, totv2)
    corr <- do_lap(R - min(R), lap_method)
    corr <- data.frame(corr_A = c(seeds$A, nonseeds$A), corr_B = c(seeds$B, nonseeds$B[corr]))
    corr <- corr[order(corr$corr_A),] 
    names(corr) <- c("corr_A","corr_B")
    rownames(corr) <- paste0(as.character(1:nrow(corr)))
    cl <- match.call()
    z <- list(call = cl, corr = corr, ns = nrow(seeds))
    z
  }
}
#' @rdname graph_match_methods
#' @return \code{graph_match_Umeyama} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B} and the number of seeds. 
#'
#' @references S. Umeyama (1988), \emph{An eigendecomposition approach to weighted
#'   graph matching problems}. IEEE TPAMI. USA, pages 695-703.
#'
#' @examples
#' # match G_1 & G_2 using Umeyama algorithm
#' G <- sample_correlated_gnp_pair(10, .9, .5)
#' G1 <- G$graph1
#' G2 <- G$graph2
#' GM_U <- graph_match_Umeyama(G1, G2, startm, alpha = .3)
#'
#' @export
#'
graph_match_Umeyama <- function(A, B, similarity = NULL, seeds = NULL, alpha = .5){

  # USE CHECK GRAPH
  A <- A[]
  B <- B[]
  totv1 <- nrow(A)
  totv2 <- nrow(B)
  
  if(totv1 > totv2){
    diff <- totv1 - totv2
    B <- pad(B[], diff)
  }else if(totv1 < totv2){
    diff <- totv2 - totv1
    A <- pad(A[], diff)
  }
  
  seeds_log <- check_seeds(seeds, nv = max(totv1, totv2), logical = TRUE)
  seeds <- check_seeds(seeds, nv = max(totv1, totv2))
  nonseeds <- seeds$nonseeds
  seeds <- seeds$seeds
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  #similarity <- similarity %*% Matrix::t(rpmat)

  if(!isSymmetric(as.matrix(A)) | !isSymmetric(as.matrix(B))){
    # construct Hermitian matrices by adjacency matrices
    A <- as.matrix((A + Matrix::t(A))/2) + as.matrix((A - Matrix::t(A))/2)*1i
    B <- as.matrix((B + Matrix::t(B))/2) + as.matrix((B - Matrix::t(B))/2)*1i
  }

  U_A <- eigen(A)$vectors
  U_B <- eigen(B)$vectors
  AB <- Matrix::tcrossprod(abs(U_B), abs(U_A))

  Grad <- alpha * AB[nonseeds$A, nonseeds$B] + (1-alpha) * Matrix::t(similarity)
  Grad <- Grad - min(Grad)
  lap_method <- set_lap_method(NULL, totv1, totv2)
  ind <- do_lap(Grad, lap_method)

  corr <- data.frame(corr_A = c(seeds$A, nonseeds$A), corr_B = c(seeds$B, nonseeds$B[ind]))
  corr <- corr[order(corr$corr_A),]
  cl <- match.call()
  z <- list(call = cl, corr = corr, ns = 0)
  z
}
