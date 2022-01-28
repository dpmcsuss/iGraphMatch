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

# sinkhorn <- function(m, niter=20){
#   # m <- matrix(abs(runif(n^2)), n)
#   for(i in 1:niter){
#     cs <- Matrix::colSums(m)
#     cs[cs == 0] <- 1
#     rs <- Matrix::rowSums(m)
#     rs[rs == 0] <- 1
#     m <- Matrix::Diagonal(x = 1 / rs) %*%
#       m %*% Matrix::Diagonal(x = 1 / cs)
#   }
#   m
# }

sinkhorn <- function(m, niter = 20) {
  d <- dim(m)
  if(d[1] != d[2]) {
    if(d[1] < d[2]) {
      return(Matrix::Diagonal(x = 1 / rowSums(m)) %*% m)
    } else {
      return(m %*% Matrix::Diagonal(x = 1 / colSums(m)))
    }
    
  }
  rs <- rep(1, d[1])
  for(i in 1:niter) {
    cs <- 1 / crossprod(m, rs)
    rs <- 1 / (m %*% cs)
  }
  Matrix::Diagonal(x = as.vector(rs)) %*% m %*%
    Matrix::Diagonal(x = as.vector(cs))
}


rds_sinkhorn <- function(n, distribution="runif"){
  sinkhorn(matrix(abs(do.call(distribution, list(n^2))), n))
}

rds_perm_bari_start <- function(nns, ns = 0, soft_seeds = NULL, g = 1, is_splr = TRUE){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  rds <- rds_perm_bari(nns - nss, g)
}



rds_perm_bari <- function(nns, g){
  alpha <- g * stats::runif(1)
  (1 - alpha) * bari_start(nns) + alpha * rperm(nns)
}

rds_from_sim_start <- function(nns, ns = 0,
    soft_seeds = NULL, sim) {

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