#' @rdname gm_fw
#' @return \code{graph_match_convex} returns a list of graph matching results,
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
  # Pseq <- list()
  while(toggle && iter < max_iter){
    f_old <- f
    iter <- iter + 1
    Grad <- -ml_sum(
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
    # Pseq <- c(Pseq, Pdir)
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
  # D <- P
  # D[nonseeds$A, nonseeds$B] <- D_ns %*% rpmat
  reorderA <- order(c(nonseeds$A, seeds$A))
  reorderB <- order(c(nonseeds$B, seeds$B))

  D <- pad(D_ns %*% rpmat, ns)[reorderA, reorderB]
  if (is(D, "splrMatrix")) {
    D@x[seeds$A, seeds$B] <- P[seeds$A, seeds$B]  
  } else {
    D[seeds$A, seeds$B] <- P[seeds$A, seeds$B]
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
  z <- list(
    call = cl,
    corr = data.frame(corr_A = 1:nv, corr_B = corr),
    ns = ns, 
    P = P,
    D = D,
    num_iter = iter
    # seq = list(alpha_seq = alpha_seq, Pseq = Pseq)
  )  
  z
}