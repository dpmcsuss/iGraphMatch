#

#' @title Solves the linear assignment problem using the Jonker-Vogenant algorithm
#'
#' @description Find a set of vertices pairs in the order of goodness of matching according to a
#' specified measure.
#'
#' @param cost A non-negative matrix-like object that can be coerced to a matrix
#' @param maximize If FALSE (default) then costs are minimized and if TRUE the
#' costs are maximized
#'
#' @details The C++ code for this method is modified from code in the
#'  \href{https://github.com/Bram94/lapjv}{python lapjv} package.
#'
#' @return The assignment of rows to columns as an integer vector
#'
lapjv <- function(cost, maximize = FALSE) {
    m <- max(cost, 10)
    n <- nrow(cost)

    cost <- rbind(cbind(as.matrix(cost), m + m * stats::runif(n)),
            m + m * stats::runif(n + 1))
    cost[n + 1, n + 1] <- 10 * m ^ 3

    ind <- cpp_lapjv(cost, maximize)
    if (ind[n + 1] <= n){
        if (sum(cost[which(ind == n + 1), 1:n]) > 1e-10){
            warning(paste("Bad padding happened. Assigned",
                which(ind == n + 1), "to", ind[n + 1]))
        }
        ind[which(ind == n + 1)] <- ind[n + 1]
    }
    ind[1:n]
}



lapmod_index <- function(n, cc, ii, kk, maximize = FALSE) {
    cpp_lapmod(n, cc, ii, kk, maximize)
}

#' @title Solves the linear assignment problem using the LAPMOD algorithm
#'
#' @description Find a set of vertices pairs in the order of goodness of matching according to a
#' specified measure.
#'
#' @param cost A non-negative CsparseMatrix object from the 'Matrix' package
#' @param maximize If FALSE (default) then costs are minimized and if TRUE the
#' costs are maximized
#'
#' @details The 'C++' code for this method is modified from code in the
#'  \href{https://github.com/Bram94/lapjv}{python lapjv} package.
#'
#' @return The assignment of rows to columns as an integer vector
#'
lapmod <- function(cost, maximize = FALSE){
    cost <- Matrix::Matrix(cost, sparse = TRUE)
    n <- nrow(cost)
    m <- max(abs(cost@x), 2)
    sign <- ifelse(maximize, -1, 1)
    pad_vec <- sign * 1e5 * m * rep(1, n)
    cost <- rbind2(cbind2(cost, sign * 1e5 * (m * ceiling(stats::runif(n))),
        c(pad_vec, - sign * 1e5 * m)))
    ind <- cpp_lapmod(n + 1, cost@x,
        cost@p, cost@i, maximize)
    if (ind[n + 1] <= n){
        if (sum(cost[which(ind == n + 1), 1:n]) > 1e-10){
            warning(paste("Bad padding happened. Assigned",
                which(ind == n + 1), "to", ind[n + 1]))
        }
        ind[which(ind == n + 1)] <- ind[n + 1]
    }
    ind[1:n]
}
