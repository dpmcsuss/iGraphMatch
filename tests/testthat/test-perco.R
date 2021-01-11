context("percolation & ExpandWhenStuck")

# sample pair of graphs w. 10 vertices
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2
seeds <- c(1, 5, 3)

test_that("perco of same sizes", {
  expect_equal(graph_match_percolation(A, B, seeds = seeds)$seeds,
               data.frame(A = c(1, 5, 3), B = c(1, 5, 3)))
})

# with similarity score
sim <- matrix(rnorm(100), 10)
test_that("perco w. similarity score", {
  expect_equal(graph_match_percolation(A, B, seeds = seeds, similarity = sim)$seeds,
               data.frame(A = c(1, 5, 3), B = c(1, 5, 3)))
})

test_that("percolation without seeds", {
  expect_equal(nrow(graph_match_percolation(A, B, seeds = NULL, similarity = sim)$seeds),
               0)
})

# directed graphs
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5, directed = TRUE)
A <- g$graph1
B <- g$graph2
test_that("perco w. similarity score", {
  expect_equal(graph_match_percolation(A, B, seeds = seeds, similarity = sim)$seeds,
               data.frame(A = c(1, 5, 3), B = c(1, 5, 3)))
})

set.seed(12)
gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
B <- lapply(gp_list, function(gp)gp[[2]])
seeds <- 1:3

test_that("percolation multi-layer", {
  expect_equal(graph_match_percolation(A, B, seeds = 1:3)$seeds,
               data.frame(A = 1:3, B = 1:3))
})



