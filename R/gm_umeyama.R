#' @title Spectral Graph Matching Methods: Umeyama Algorithm
#' @rdname gm_Umeyama
#' 
#' @param A A matrix, igraph object, or list of either.
#' @param B A matrix, igraph object, or list of either. 
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similaities.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix.
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#'   
#' @return \code{graph_match_Umeyama} returns a list of graph matching 
#'   results, including the graph matching formula, a data frame containing the 
#'   matching correspondence between \eqn{G_1} and \eqn{G_2} named \code{corr_A} 
#'   and \code{corr_B} and seeds. 
#'
#' @references S. Umeyama (1988), \emph{An eigendecomposition approach to weighted
#'   graph matching problems}. IEEE TPAMI. USA, pages 695-703.
#'
#' @examples
#' # match G_1 & G_2 using Umeyama algorithm
#' G <- sample_correlated_gnp_pair(10, .9, .5)
#' g1 <- G$graph1
#' g2 <- G$graph2
#' startm <- matrix(0, 10, 10)
#' diag(startm)[1:4] <- 1
#' graph_match_Umeyama(g1, g2, similarity = startm)
#'
#' @export
#'
graph_match_Umeyama <- function(A, B, seeds = NULL, 
                                similarity = NULL){
  
  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  nc <- length(A)
  
  seeds <- check_seeds(seeds, nv = max(totv1, totv2))
  nonseeds <- seeds$nonseeds
  seeds <- seeds$seeds
  similarity <- check_sim(similarity, seeds, nonseeds, totv1, totv2)
  Grad <- 0
  
  for( ch in 1:nc ){
    if(!isSymmetric(as.matrix(A[[ch]])) | !isSymmetric(as.matrix(B[[ch]]))){
      # construct Hermitian matrices by adjacency matrices
      A[[ch]] <- as.matrix((A[[ch]] + Matrix::t(A[[ch]]))/2) + as.matrix((A[[ch]] - Matrix::t(A[[ch]]))/2)*1i
      B[[ch]] <- as.matrix((B[[ch]] + Matrix::t(B[[ch]]))/2) + as.matrix((B[[ch]] - Matrix::t(B[[ch]]))/2)*1i
    }
    
    U_A <- eigen(A[[ch]])$vectors
    U_B <- eigen(B[[ch]])$vectors
    AB <- Matrix::tcrossprod(abs(U_B), abs(U_A))
    Grad <- Grad + AB[nonseeds$A, nonseeds$B]
  }
  

  Grad <- Grad + Matrix::t(similarity)
  Grad <- Grad - min(Grad)
  
  # make a random permutation
  nn <- nrow(A[[1]]) - nrow(seeds)
  rp <- sample(nn)
  rpmat <- Matrix::Diagonal(nn)[rp, ]
  Grad <- Grad[,rp]
  
  lap_method <- set_lap_method(NULL, totv1, totv2)
  ind <- do_lap(Grad, lap_method)
  
  # undo rand perm here
  ind <- rp[ind]
  corr <- data.frame(corr_A = c(seeds$A, nonseeds$A), corr_B = c(seeds$B, nonseeds$B[ind]))
  corr <- corr[order(corr$corr_A),]
  cl <- match.call()
  z <- list(
    call = cl, 
    corr = corr, 
    seeds = seeds)
  z
}
