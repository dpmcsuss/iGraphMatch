#' @title Start matrix initialization
#'
#' @description initialize the start matrix for graph matching iteration.
#'
#' @param nns An integer. Number of non-seeds.
#' @param ns An integer. Number of hard seeds.
#' @param soft_seeds A vector, a matrix or a data frame. If there is no error in soft
#' seeds, input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft
#' seeds, input in the form of a matrix or a data frame, with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}. Note
#' that if there are seeds in graphs, seeds should be put before non-seeds.
#' @param distribution A charactor. Specify the distribution from which the random doubly stochastic
#' matrix is sampled. Should input the name of the function for generating random deviates from that
#' distribution.
#' @param g A number. Specified in the range of [0, 1] to set weights to random permutaion matrix and
#' barycenter matrix.
#'
#' @rdname start
#' @return \code{bari_start} returns a \code{nns-by-nns} matrix with 1's corresponding to the
#' adaptive seeds and being bari-centered at other places.
#' @examples
#' ## Case without soft seeds
#' bari_start(3)
#'
#' ## Case with correct soft seeds and input is a vector
#' bari_start(nns=5, ns=3, soft_seeds=c(5, 7, 8))
#'
#' ## Case with erroneous soft seeds and the input is a matrix
#' bari_start(nns=5, soft_seeds=matrix(c(2, 4, 2, 3), nrow=2))
#'
#' @export
#'
bari_start <- function(nns, ns = 0, soft_seeds = NULL){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  start <- bari_splr(nns - nss)
  add_soft_seeds(start, nns, ns, soft_seeds)
}

bari_splr <- function(nns){
  nns <- as.integer(nns)
  splr(x = Matrix(0, nns, nns), 
    a = Matrix(1, nns), b = Matrix(1 / nns, nns))
}

#' @rdname start
#' @return \code{rds_sinkhorn_start} returns a \code{nns-by-nns} doubly stochastic matrix
#' with 1's corresponding to adaptive seeds.
#' @examples
#' ## Case without soft seeds
#' rds_sinkhorn_start(5)
#'
#' ## Case with soft seeds and the input is a data frame
#' rds_sinkhorn_start(nns=5,
#'    soft_seeds = as.data.frame(matrix(c(2, 4, 2, 3), nrow=2)),
#'    distribution = "rnorm")
#' 
#' @export
#'
rds_sinkhorn_start <- function(nns, ns = 0, soft_seeds = NULL, distribution = "runif"){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  rds <- rds_sinkhorn(nns - nss,
    distribution = distribution)
  add_soft_seeds(rds, nns, ns, soft_seeds)
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

 
#' @rdname start
#' 
#' @param is_splr should we return a splr matrix? (default = TRUE)
#' 
#' @return \code{rds_perm_bari} returns a \code{nns-by-nns} doubly stochastic matrix
#' with 1's corresponding to adaptive seeds.
#' 
#' @examples
#' ## Case without soft seeds
#' rds_perm_bari_start(nns=5)
#'
#' ## Case with soft seeds and the input is a data frame
#' rds_perm_bari_start(nns=5, ns=0, soft_seeds=as.data.frame(matrix(c(2, 4, 2, 3), nrow=2)))
#' 
#' @export
#' 
rds_perm_bari_start <- function(nns, ns = 0, soft_seeds = NULL, g = 1, is_splr = TRUE){
  nss <- nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  rds <- rds_perm_bari(nns - nss, g, is_splr)
  add_soft_seeds(rds, nns, ns, soft_seeds)
}



rds_perm_bari <- function(nns, g, is_splr = TRUE){
  alpha <- stats::runif(1, 0, g)
  if(is_splr){
    if(is.numeric(alpha * rperm(nns))){
      x <- sparseMatrix(x = alpha * rperm(nns), i=1, j=1)
    } else{
      x <- alpha * rperm(nns)
    }
    new("splrMatrix",
        x = x, 
        a = Matrix(1 - alpha, nns), b = Matrix(1 / nns, nns),
        Dim = c(as.integer(nns), as.integer(nns)),
        Dimnames = list(NULL, NULL))
  } else {
    (1 - alpha) * bari_start(nns) +
        alpha * rperm(nns)
  }
}

#' @rdname start
#' 
#' @param sim nns x nns non-negative matrix.
#' 
#'
#' @return \code{rds_from_sim_start} returns a doubly
#'  stochastic Matrix given by sinkhorn algorithm applied to 
#'  a matrix of iid log-normal with mu=sim. Note,
#'  this ignores soft seeds.
#' 
#' @examples
#' sim <- Matrix::rsparsematrix(10, 10, .4,
#'  rand.x = function(n) rep(1,n))
#' start_sparse <- rds_from_sim_start(10, sim)
#' start_dense <- rds_from_sim_start(10, as.matrix(sim))
#' 
#' @export
rds_from_sim_start <- function(nns, ns = 0,
    soft_seeds = NULL, sim) {

  if (!is.null(soft_seeds)) {
    warning("Ignoring soft_seeds in rds_from_sim_start")
  }
  rds_from_sim(nns, sim)
}

rds_from_sim <- function(nns, sim) {
  if (inherits(sim, "sparseMatrix") &&
      "x" %in% slotNames(sim)) {
    sim@x <- exp(sim@x + stats::rnorm(Matrix::nnzero(sim)) * 2)
    sinkhorn(sim, 40)
  } else {
    sinkhorn(Matrix(exp(stats::rnorm(nns ^ 2, 1)), nns) + sim, 40)
  }
}



# #' @rdname featurestart
# #' 
# #' @param f1
# #' @param f2
# #' @param s
# #' 
# #' @return \code{feature_tart} returns a 
# #' \code{nns-by-nns} matrix
# #' with entries corresponding to similiarity scores
# #' for feature vectors.
# #' 
# #' @examples
# feature_start <- function(f1, f2, s){

# }