#' @title Select adaptive seeds
#'
#' @description Find a set of vertex pairs that are best matched be to adaptive seeds.
#'
#' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
#' @param measure A character. Measure for computing goodness of matching.
#' @param nseeds An integer. Number of adaptive seeds needed.
#' @param x A vector of logical. \code{TRUE} indicates the corresponding vertex is of interest
#' in finding adaptive seeds. Length of vector should be the number of vertice of graphs.
#' @param match_corr A vector. Graph matching correspondence between \eqn{G_1} and \eqn{G_2}.
#'
#' @return \code{select_seeds} returns a data frame with the indice of seeds in \eqn{G_1} named
#' \code{seed_A} and the indice of seeds in \eqn{G_2} named \code{seed_B}.
#'
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 50, corr =  0.3, p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' seeds <- 1:50 <= 10
#' nonseeds <- !seeds
#' match <- graph_match_FW(g1, g2, seeds)
#'
#' # select adaptive seeds from non seeds
#' select_seeds(g1, g2, "row_perm_stat", nseeds = 5, x = nonseeds, match$corr)
#'
#' @export
#'
select_seeds <- function(A, B, measure, nseeds, x, match_corr){
  A <- A[]
  B <- B[]
  A <- as.matrix(A)
  B <- as.matrix(B)
  Bm <- B[match_corr, match_corr]
  nv <- dim(A)[1]

  # calculate measure stat
  stat <- do.call(measure,list(A,Bm))
  stat_nonseeds <- stat[x]

  # find top ranking nodes pairs
  rstat <- sort(stat_nonseeds)
  topstat <- rstat[1:nseeds]
  topindex <- sapply(topstat, function(top) which(stat==top)) %>%
    unlist %>% unique # solve ties in topstat
  topindex <- sapply(topindex, function(top){ifelse(top %in% c(1:nv)[!x],0,top)})
  topindex <- topindex[topindex!=0]
  topindex <- topindex[1:nseeds]

  match_nseeds <- match_corr[topindex]

  selected_seeds <- data_frame(seed_A=topindex, seed_B=match_nseeds)
  selected_seeds
}
