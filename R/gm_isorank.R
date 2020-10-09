#' @title Spectral Graph Matching Methods: IsoRank Algorithm
#' @rdname gm_isorank
#' 
#' @param A A matrix, igraph object, or list of either.
#' @param B A matrix, igraph object, or list of either. 
#' @param similarity A matrix. An \code{n-by-n} matrix containing vertex similaities.
#' @param seeds A vector of integers or logicals, a matrix or a data frame. If
#'   the seed pairs have the same indices in both graphs then seeds can be a
#'   vector. If not, seeds must be  a matrix
#'   or a data frame, with the first column being the indices of \eqn{G_1} and
#'   the second column being the corresponding indices of \eqn{G_2}.
#' @param alpha A number betwen 0 and 1. Bigger alpha means putting more importance
#'   on the information in network topology over other information such as
#'   similarity scores
#' @param max_iter A number. Maximum number of replacing matches equals to
#'   max_iter times number of total vertices of \eqn{G_1}.
#' @param method A character. Choice of method to extract mapping from score matrix,
#'   including greedy method and the Hungarian algorithm.
#' 
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
#' cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' # match G_1 & G_2 using IsoRank algorithm
#' startm <- matrix(0, 10, 10)
#' diag(startm)[1:4] <- 1
#' GM_IsoRank <- graph_match_IsoRank(g1, g2, startm, alpha = .3, method = "greedy")
#'
#' @export
#'
graph_match_IsoRank <- function(A, B, similarity, seeds = NULL, 
                                alpha = .5, max_iter = 50, method = "greedy"){
  
  graph_pair <- check_graph(A, B)
  A <- graph_pair[[1]]
  B <- graph_pair[[2]]
  totv1 <- graph_pair$totv1
  totv2 <- graph_pair$totv2
  nc <- length(A)
  
  seeds <- check_seeds(seeds, nv = max(totv1, totv2))
  nonseeds <- seeds$nonseeds
  seeds <- seeds$seeds
  if(!is.null(similarity) && dim(similarity)[1] != dim(similarity)[2]){
    diff <- dim(similarity)[1] - dim(similarity)[2]
    similarity <- pad(similarity, max(-diff, 0), max(diff, 0))
  }
  R <- E <- similarity / sum(abs(similarity))
  tol <- 1e-2
  R_tot <- matrix(0, nrow(R), ncol(R))
  
  for( ch in 1:nc ){
    
    iter <- 1
    diff <- 1
    
    # computing transition matrix A
    colS_A <- Matrix::colSums(A[[ch]])
    colS_B <- Matrix::colSums(B[[ch]])
    A[[ch]] <- A[[ch]] %*% Matrix::Diagonal(nrow(A[[ch]]), ifelse(colS_A == 0, 0, 1/colS_A))
    B[[ch]] <- B[[ch]] %*% Matrix::Diagonal(nrow(B[[ch]]), ifelse(colS_B == 0, 0, 1/colS_B))
    
    # computing R by power method
    while(diff > tol & iter <= max_iter){
      
      if(alpha>0){
        AR <- A[[ch]] %*% R %*% Matrix::t(B[[ch]])
        AR <- alpha * AR + (1-alpha) * E
      } else{
        AR <- A[[ch]] %*% R %*% Matrix::t(B[[ch]])
      }
      R_new <- AR / sum(abs(AR))
      diff <- sum(abs(R-R_new))
      iter <- iter + 1
      R <- R_new
    }
    
    R <- as.matrix(R)
    R_tot <- R_tot + R
  }
  
  # find GNA
  R <- R_tot[nonseeds$A, nonseeds$B]
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