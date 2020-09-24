#' @rdname gm_spectral
#' 
#' @title Spectral Graph Matching Methods
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
