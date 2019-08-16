context("percolation & ExpandWhenStuck")

# sample pair of graphs w. 10 vertices
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2
seeds <- c(1, 5, 3)

test_that("perco & ExpandWhenStuck of same sizes", {
  expect_equal(graph_match_percolation(A, B, seeds = seeds)$ns, 3)
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds)$ns, 3)
})

