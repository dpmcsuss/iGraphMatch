#

#' @title Solves a linear assignment problem using the Jonker-Vogenant algorithm or LAPMOD variant
#'
#' @description Find the matching of rows to columns that minimizes or maximizes 
#' the cost. See \link{do_lap} for usage.
#'
#' @param cost For lapjv, an object that can be coerced to a matrix. For lapmod, a sparseMatrix.
#' @param maximize If FALSE (default) then costs are minimized and if TRUE the
#' costs are maximized
#'
#' @details The C++ code for these method is modified from code in the
#'  \href{https://github.com/Bram94/lapjv}{python lapjv} package.
#'  
#'  The cost matrix is padded with a single row and column of very large entries that
#'  helps to avoid stability issues with the algorithms.
#'
#' @return The assignment of rows to columns as a vector.
#'
#' @rdname lapjv
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

#' @rdname lapjv
lapmod <- function(cost, maximize = FALSE){
    cost <- Matrix::Matrix(cost, sparse = TRUE)
    n <- nrow(cost)
    m <- max(abs(cost@x), 2)
    sign <- ifelse(maximize, -1, 1)
    pad_vec <- sign * 1e5 * m * rep(1, n)
    cost <- 
        rbind(
            cbind(
                cost,
                sign * ceiling(m * 1e5 * stats::runif(n))
            ),
            c(pad_vec, - sign * 1e5 * m)
        )
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
