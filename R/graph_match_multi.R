
#' @title Multiple Graph Match Methods
#'
#' @description Match two lists of graphs, returns a list of graph matching results,
#' including matching correspondence vector of \eqn{G_2} with respect to \eqn{G_1},
#' doubly stochastic matrix and permutation matrix.
#'
#' @param A A list of graphs (or igraph objects)
#' @param B A list of graphs (or igraph objects)
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If there is no error in seeds input can be
#' a vector of seed indices in \eqn{G_1}. Or if there exists error in seeds, input in the form of a
#' matrix or a data frame, with the first column being the indices of \eqn{G_1} and the second
#' column being the corresponding indices of \eqn{G_2}.
#' @param start A matrix or a character. Any \code{nns-by-nns} matrix or
#' character value like "bari" or "convex" to initialize the starting matrix.
#' @param max_iter An integer. Maximum iteration time.
#' @param tol A number. Tolerance of edge disagreements.
#' @param r A number. Threshold of neighboring pair scores.
#'
#' @rdname graph_match_methods
#'
#' @return \code{graph_match_FW} returns a list of graph matching results,
#' including match correspondence vector of \eqn{G_2} with respect to \eqn{G_1}
#' named \code{corr}, doubly stochastic matrix named \code{D}, permutation
#' matrix named \code{P} based on Frank-Wolfe methodology and iteration time of
#' the algorithm named \code{iter}.
#'
#' @examples
#'  gp_list <- replicate(3, sample_correlated_gnp_pair(100, .3, .5), simplify = FALSE)
#'  A <- lapply(gp_list, function(gp)gp[[1]])
#'  B <- lapply(gp_list, function(gp)gp[[2]])
#'  match <- graph_match_FW_multi(A, B, seeds = 1:10, start = "bari", max_iter = 20)
#'  match$corr
#'
#' @export
#'
graph_match_FW_multi <- function(A, B, seeds = NULL, start = "bari", max_iter = 20){

  if(start == "convex"){
    stop("Convex start is not yet implemented for multiplex matching")
  }

  # this will make the graphs be matrices if they are igraph objects
  A <- lapply(A, function(Al) Al[])
  B <- lapply(B, function(Bl) Bl[])

  # Assume each list has all the same number of nodes within
  totv1<-ncol(A[[1]])
  totv2<-ncol(B[[1]])
  if(totv1>totv2){
    diff<-totv1-totv2
    B <- lapply(B, function(Bl)
      Matrix::bdiag(Bl[], Matrix(0,diff,diff)))
  }else if(totv1<totv2){
    diff<-totv2-totv1
    A <- lapply(A, function(Al)
      Matrix::bdiag(Bl[], Matrix(0,diff,diff)))
  }
  nv <- nrow(A[[1]])

  if(is.null(seeds)){
    seeds <- rep(FALSE,nv)
    aseeds_err <- FALSE
    ns <- sum(seeds)
  } else{
    seeds_pair <- check_seeds(seeds)
    ns <- nrow(seeds_pair)

    seeds <- rep(FALSE,nv)
    seeds[seeds_pair$seed_A] <- TRUE

    # detect incorrect seeds
    seed_A <- seeds_pair$seed_A
    seed_B <- seeds_pair$seed_B
    aseeds_err <- ifelse(seed_A!=seed_B,TRUE,FALSE)
    seed_A_err <- seed_A[aseeds_err]
    seed_B_err <- seed_B[aseeds_err]

    if(sum(aseeds_err)!=0){
      B <- g2_hard_seeding(seed_A_err,seed_B_err,B)
    }
  }

  nn <- nv-ns
  nonseeds <- !seeds


  P <- init_start(start = start, nns = nn,
    A = A, B = B, seeds = seeds)

  iter <- 0
  toggle <- TRUE

  # seed to non-seed info
  s_to_ns <- get_s_to_ns(A,B, seeds)

  # keep only nonseeds
  A <- lapply(A, function(Al) Al[nonseeds, nonseeds])
  B <- lapply(B, function(Bl) Bl[nonseeds, nonseeds])

  while(toggle && iter < max_iter){
    iter <- iter + 1

    # non-seed to non-seed info
    tAnn_P_Bnn <- Reduce("+", mapply(
      function(Ann, Bnn){ Matrix::t(Ann) %*% P %*% Bnn},
      A, B))

    Grad <- s_to_ns + tAnn_P_Bnn + Reduce("+", mapply(
        function(Ann, Bnn){ Ann %*% P %*% Matrix::t(Bnn)},
      A, B))
    Grad <- Grad-min(Grad)

    ind <- as.vector(clue::solve_LSAP(as.matrix(Grad), maximum = TRUE))
    ind2 <- cbind(1:nn, ind)
    Pdir <- Matrix::Diagonal(nn)
    Pdir <- Pdir[ind, ]
    ns_Pdir_ns <- Reduce("+", mapply(
        function(Ann, Bnn) Matrix::t(Ann)[, order(ind)] %*% Bnn,
        A, B))
    c <- sum(tAnn_P_Bnn * P)
    d <- sum(ns_Pdir_ns * P) + sum(tAnn_P_Bnn[ind2])
    e <- sum(ns_Pdir_ns[ind2])
    u <- sum(P * (s_to_ns))
    v <- sum((s_to_ns)[ind2])
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
  corr_ns <- as.vector(clue::solve_LSAP(as.matrix(round(P*nn^2)), maximum = TRUE))
  corr <- 1:nv
  corr[nonseeds] <- corr[nonseeds][corr_ns]
  P <- Matrix::Diagonal(nv)[corr,]
  D <- P
  D[nonseeds,nonseeds] <- D_ns

  # fix match results if there are incorrect seeds
  if(sum(aseeds_err)!=0){
    corr <- fix_hard_corr(seed_A_err,seed_B_err,corr)
    P <- Matrix::Diagonal(nv)[corr,]
    D <- fix_hard_D(seed_A_err,seed_B_err,D)
  }

  list(corr = corr, P = P, D = D, iter = iter)
}

get_s_to_ns <- function(Alist, Blist, seeds){
  nonseeds <- !seeds
  s_to_ns <- function(A,B){
    Asn <- A[seeds,nonseeds]
    Ann <- A[nonseeds,nonseeds]
    Ans <- A[nonseeds,seeds]

    Bsn <- B[seeds,nonseeds]
    Bnn <- B[nonseeds,nonseeds]
    Bns <- B[nonseeds,seeds]
    Ans %*% Matrix::t(Bns) + Matrix::t(Asn) %*% Bsn
  }
  Reduce("+", mapply(s_to_ns, Alist, Blist, SIMPLIFY = FALSE))
}
