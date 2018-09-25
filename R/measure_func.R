#' @title Measure functions
#'
#' @description Measures for computing the goodness of matching for each vertex.
#'
#' @param g1 A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
#' @param g2 A matrix or an igraph object. Adjacency matrix of \eqn{G_2} after adjusting rows and columns
#' according to the correlation of matching between two graphs.
#' @param exact A logical. If \code{g1} and \code{g2} are binary, then set \code{exact=TRUE},
#' if \code{g1} and \code{g2} are weighted graphs, then set \code{exact=FALSE}.
#'
#'
#' @rdname measure_func
#' @return \code{row_cor} returns a vector, each element is 1 minus the row correlation value for
#' the corresponding vertex.
#' @examples
#' cgnp_pair <- sample_correlated_gnp_pair(n = 50, corr =  0.3,p =  0.5)
#' g1 <- cgnp_pair$graph1
#' g2 <- cgnp_pair$graph2
#' match <- graph_match_FW(g1, g2)
#' g2m <- g2[match$corr, match$corr]
#' g1 <- g1[]
#' row_cor(g1, g2m)
#' @export
#'
row_cor <- function(g1,g2){
  g1 <- g1[]
  g2 <- g2[]

  require(tidyverse)
  1:nrow(g1) %>% map_dbl(~suppressWarnings(1-cor(g1[.x,],g2[.x,])))
}
#'
#' @rdname measure_func
#' @return \code{row_diff} returns a vector, each element is the row difference value for
#' the corresponding vertex.
#' @examples
#' row_diff(g1, g2m)
#' @export
#'
row_diff <- function(g1,g2){
  g1 <- g1[]
  g2 <- g2[]
  g1 <- as.matrix(g1)
  g2 <- as.matrix(g2)
  rowSums(abs(g1-g2))
}
#'
#' @rdname measure_func
#' @return \code{row_perm_stat} returns a vector, each element is the row permutation statistics
#' value for the corresponding vertex.
#' @examples
#' row_perm_stat(g1, g2m)
#' @export
#'
row_perm_stat <- function(g1,g2,exact=TRUE,...){
  g1 <- g1[]
  g2 <- g2[]

  if(exact){
    m <- mean_row_diff(g1,g2)
    v <- var_row_diff(g1,g2)
  } else {
    mv <- row_diff_perm(g1,g2)
    m <- mv$mean
    v <- mv$var
  }

  d <- rowSums(abs(g1-g2))

  (d-m)/sqrt(v)
}

row_diff_perm <- function(g1, g2, nmc = 1000, sym=FALSE){
  g1 <- g1[]
  g2 <- g2[]

  n <- nrow(g2)
  A <- Matrix(0,n,nmc)

  for(mc in 1:nmc){
    p <- sample(n)
    A[,mc] <- rowSums(abs(g1[,order(p)]-g2[p,]))
  }
  m <- rowMeans(A)
  v <- apply(A,1,var)
  if(sym){
    mv <- row_diff_perm(g1,g2,nmc)
    m <- m+mv$mean
    v <- sqrt(v*mv$var)
  }
  list(mean=m,var=v)
}
mean_row_diff <- function(g1, g2, sym=FALSE){
  g1 <- g1[]
  g2 <- g2[]

  dg1 <- rowSums(g1)
  dg2 <- rowSums(g2)
  mdg2 <- mean(dg2)
  n <- nrow(g1)

  r1 <- mdg2
  r2 <- dg1 * (1 - 2 * mdg2 / (n - 1))
  ED <- r1 + r2
  if (sym){
    ED <- (ED + mean_row_diff(g2, g1)) / 2
  }
  ED
}
var_row_diff <- function(g1, g2, sym=FALSE){
  g1 <- g1[]
  g2 <- g2[]

  dg1 <- rowSums(g1)
  dg2 <- rowSums(g2)
  mdg2 <- mean(dg2)
  n <- nrow(g1)

  m2dg2 <- mean(dg2^2)
  vdg2 <- m2dg2-mdg2^2 # don't use var here

  v1 <- (1-2*dg1/(n-1))^2*vdg2
  v2 <- 4*dg1*(n-1-dg1)*((n-1)*mdg2-m2dg2)/((n-1)^2*(n-2))
  VD <- v1+v2

  if(sym){
    VD <- sqrt(VD*var_row_diff(g2,g1))
  }
  VD
}
