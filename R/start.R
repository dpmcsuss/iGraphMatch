bari_start <- function(nns, ns = 0, soft_seeds = NULL){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  start <- bari_splr(nns - nss)
}

bari_splr <- function(nns){
  nns <- as.integer(nns)
  splr(x = Matrix(0, nns, nns),
    a = Matrix(1, nns), b = Matrix(1 / nns, nns))
}

rds_sinkhorn_start <- function(nns, ns = 0, soft_seeds = NULL, distribution = "runif", ...){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  rds <- rds_sinkhorn(nns - nss,
    distribution = distribution)
}

sinkhorn <- function(m, niter=20){
  # m <- matrix(abs(runif(n^2)), n)
  for(i in 1:niter){
    r <- Matrix::colSums(m)
    r[r == 0] <- 1
    m <- t(m %*% Matrix::Diagonal(x = 1 / r))
  }
  m
}

rds_sinkhorn <- function(n, distribution="runif"){
  sinkhorn(matrix(abs(do.call(distribution, list(n^2))), n))
}

rds_perm_bari_start <- function(nns, ns = 0, soft_seeds = NULL, g = 1, is_splr = TRUE){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  rds <- rds_perm_bari(nns - nss, g)
}



rds_perm_bari <- function(nns, g, ...){
  alpha <- g * stats::runif(1)
  (1 - alpha) * bari_start(nns) + alpha * rperm(nns)
}

rds_from_sim_start <- function(nns, ns = 0,
    soft_seeds = NULL, sim, ...) {

  if (!is.null(soft_seeds) && nrow(soft_seeds) > 0) {
    warning("Ignoring soft_seeds in rds_from_sim_start")
  }
  if (!is(sim, "matrix") && !is(sim, "Matrix")) {
    stop(
      paste0(
        "Error: sim must be a matrix-like object, not a ",
        class(sim)
      )
    )
  }
  rds_from_sim(nns, sim)
}

rds_from_sim <- function(nns, sim) {
  if (inherits(sim, "sparseMatrix") &&
      "x" %in% slotNames(sim)) {
    sim@x <- exp(sim@x + stats::rnorm(Matrix::nnzero(sim)) * 2)
    sim <- sinkhorn(sim, 40)

  } else {
    sim <- sinkhorn(
      Matrix(exp(stats::rnorm(prod(dim(sim)), 1)), nns) + sim,
      40
    )
  }
  diff <- max(dim(sim)) - dim(sim)
  sim <- pad(sim, diff[1], diff[2])
  sim
}