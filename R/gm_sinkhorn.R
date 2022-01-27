#' @rdname gm_sinkhorn
#'
#' @references Y. Aflalo and A. Bronstein and R. Kimmel (2014), \emph{On convex
#' relaxation of graph isomorphism}. Proceedings of the National Academy of Sciences,
#' pages 2942-2947.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 500, corr =  0.9, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 with no seeds
#' gm(g1, g2, method = "sinkhorn", max_iter = 10)
#' seeds <- 1:10 <= 3
#' gm(g1, g2, seeds, method = "sinkhorn", max_iter = 10)
#'
#'
#' @keywords internal
graph_match_sinkhorn <- function(A, B, seeds = NULL,
  similarity = NULL, start = "bari", max_iter = 100,
  lambda = 100,
  tol = 1e-5) {

  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  nv <- max(totv1, totv2)
  nonseeds <- check_seeds(seeds, nv)$nonseeds
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

  alpha_seq <- NULL
  # Pseq <- list()
  while(toggle && iter < max_iter){
    f_old <- f
    iter <- iter + 1
    Grad <- -as.matrix(ml_sum(
      AtA %*% P + P %*% BBt - ABns_sn +
      -t(Ann) %*% P %*% Bnn - Ann %*% P %*% t(Bnn) +
      similarity))

    Pdir <- sinkhorn(exp(lambda * Grad / max(Grad)))


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

    aopt <- ifelse(aq == 0 & bq == 0, 0,
      ifelse(-bq / aq > 1, 1, -bq / aq))
    alpha_seq <- c(alpha_seq, aopt)
    # Pseq <- c(Pseq, Pdir)
    P_new <- aopt * P + (1 - aopt) * Pdir
    f <- innerproduct(Ann %*% P_new - P_new %*% Bnn,
      Ann %*% P_new - P_new %*% Bnn)

    f_diff <- abs(f - f_old)
    P_diff <- norm(P - P_new, "f")
    P <- P_new
    
    toggle <- f_diff > tol && f > tol && P_diff > tol
  }

  if(iter == max_iter){
    warning("Frank-Wolfe iterations reach the maximum iteration, convergence may not occur.")
  }


  corr_ns <- do_lap(P, "lapjv") # change to something better
  # undo rand perm here
  corr_ns <- rp[corr_ns]

  corr <- 1:nv
  corr[nonseeds$A] <- nonseeds$B[corr_ns]
  corr[seeds$A] <- seeds$B
  # P <- Matrix::Diagonal(nv)[corr, ]

  reorderA <- order(c(nonseeds$A, seeds$A))
  reorderB <- order(c(nonseeds$B, seeds$B))

  D <- pad(P %*% rpmat, ns)[reorderA, reorderB]
  if (is(D, "splrMatrix")) {
    D@x[seeds$A, seeds$B] <- Matrix::Diagonal(ns)
      # P[seeds$A, seeds$B]
  } else {
    D[seeds$A, seeds$B] <- Matrix::Diagonal(ns)
     # P[seeds$A, seeds$B]
  }


  # get_f <- function(a){
  #   PP <- a * P + (1 - a) * Pdir
  #   innerproduct(Ann %*% PP - PP %*% Bnn,
  #     Ann %*% PP - PP %*% Bnn)
  # }
  # tibble(a = seq(0, 1, 0.01)) %>%
  #   mutate(f = map_dbl(a, get_f)) %>%
  #   ggplot(aes(a,f)) + geom_point() +
  #   geom_vline(xintercept = aopt)

  cl <- match.call()

  graphMatch(
    corr = data.frame(corr_A = 1:nv, corr_B = corr),
    nnodes = c(totv1, totv2),
    call = cl,
    detail = list(
      iter = iter,
      max_iter = max_iter,
      seeds = seeds,
      soft = D
    )
  )
}
