#' @title Start matrix initialization
#'
#' @description initialize the start matrix for graph matching iteration.
#'
#' @param nns An integer. Number of non-seeds.
#' @param ns An integer. Number of seeds.
#' @param soft_seeds A vector, a matrix or a data frame. If there is no error in soft
#' seeds, input can be a vector of soft seed indices in \eqn{G_1}. Or if there is error in soft
#' seeds, input in the form of a matrix or a data frame, with the first column being the
#' indices of \eqn{G_1} and the second column being the corresponding indices of \eqn{G_2}. Note
#' that if there are seeds in graphs, seeds should be put before non-seeds.
#' @param distribution A charactor. Specify the distribution from which the random doubly stochastic
#' matrix is sampled. Should input the name of the function for generating random deviates from that
#' distribution.
#' @param g A number. Specified in the range of [0,1] to set weights to random permutaion matrix and
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
#' bari_start(nns=5,ns=3,soft_seeds=c(5,7,8))
#'
#' ## Case with erroneous soft seeds and the input is a matrix
#' bari_start(nns=5,soft_seeds=matrix(c(2,4,2,3),nrow=2))
#'
#' @export
#'
bari_start <- function(nns, ns = 0, soft_seeds = NULL){
  if(is.null(soft_seeds)){
    start <- matrix(1/nns,nns,nns)
  } else{
    soft_seeds <- check_seeds(soft_seeds)
    seed_g1 <- soft_seeds$seed_A
    seed_g2 <- soft_seeds$seed_B
    nseeds <- length(seed_g1)

    start <- matrix(1/(nns-nseeds),nns,nns)
    for (i in 1:nseeds) {
      start[seed_g1[i]-ns,] <- 0
      start[,seed_g2[i]-ns] <- 0
      start[seed_g1[i]-ns,seed_g2[i]-ns] <- 1
    }
  }

  start
}

#' @rdname start
#' @return \code{rds_sinkhorn_start} returns a \code{nns-by-nns} doubly stochastic matrix
#' with 1's corresponding to adaptive seeds.
#' @examples
#' ## Case without soft seeds
#' rds_sinkhorn_start(5)
#'
#' ## Case with soft seeds and the input is a data frame
#' rds_sinkhorn_start(nns=5, soft_seeds=as.data.frame(matrix(c(2,4,2,3),nrow=2)), distribution = "rnorm")
#' @export
#'
rds_sinkhorn_start <- function(nns, ns = 0, soft_seeds = NULL, distribution = "runif"){
  if(is.null(soft_seeds)){
    start <- rds_sinkhorn(nns,distribution = distribution)
  } else{
    soft_seeds <- check_seeds(soft_seeds)
    seed_g1 <- soft_seeds$seed_A
    seed_g2 <- soft_seeds$seed_B
    nseeds <- length(seed_g1)
    
    start <- matrix(5,nrow = nns,ncol = nns)
    for (i in 1:nseeds) {
      start[seed_g1[i]-ns,] <- 0
      start[,seed_g2[i]-ns] <- 0
      start[seed_g1[i]-ns,seed_g2[i]-ns] <- 1
    }

    rds <- rds_sinkhorn(nns-nseeds,distribution = distribution)
    start[start==5] <- rds
  }

  start
}

sinkhorn <- function(m,niter=20){
  # m <- matrix(abs(runif(n^2)),n)
  for(i in 1:niter){
    r <- rowSums(m)
    m <- t(diag(1/r) %*% m)
  }
  m
}
rds_sinkhorn <- function(n,distribution="runif"){
  sinkhorn(matrix(abs(do.call(distribution,list(n^2))),n))
}
#' @rdname start
#' @return \code{rds_perm_bari} returns a \code{nns-by-nns} doubly stochastic matrix
#' with 1's corresponding to adaptive seeds.
#' @examples
#' ## Case without soft seeds
#' rds_perm_bari_start(nns=5)
#'
#' ## Case with soft seeds and the input is a data frame
#' rds_perm_bari_start(nns=5, ns=0, soft_seeds=as.data.frame(matrix(c(2,4,2,3),nrow=2)))
#' @export
#' 
rds_perm_bari_start <- function(nns, ns = 0, soft_seeds = NULL, g = 1){
  if(is.null(soft_seeds)){
    start <- rds_perm_bari(nns, g = g)
  } else{
    soft_seeds <- check_seeds(soft_seeds)
    seed_g1 <- soft_seeds$seed_A
    seed_g2 <- soft_seeds$seed_B
    nseeds <- length(seed_g1)
    
    start <- matrix(5,nrow = nns,ncol = nns)
    for (i in 1:nseeds) {
      start[seed_g1[i]-ns,] <- 0
      start[,seed_g2[i]-ns] <- 0
      start[seed_g1[i]-ns,seed_g2[i]-ns] <- 1
    }

    rds <- rds_perm_bari(nns-nseeds, g = g)
    start[start==5] <- rds
  }

  start
}


rds_perm_bari <- function(nns, g){
    alpha <- runif(1, 0, g)
    (1 - alpha) * bari_start(nns) +
        alpha * rperm(nns)
}
