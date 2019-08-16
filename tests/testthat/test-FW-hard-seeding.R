context("hard-seeding Frank Wolfe with incorrect seeds")

# sample pair of graphs w. 10 vertices
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2

# sample hard seeds
hard_cor <- c(1, 5, 3) # correct
hard_incor <- data.frame(A = c(1, 2, 3), B = c(1, 2, 5)) # contain incorrect

test_that("FW with no seed", {
  expect_equal(graph_match_FW(A, B, seeds = NULL, start = "bari")$ns, 0)
})

test_that("FW with correct hard seeds", {
  expect_equal(graph_match_FW(A, B, seeds = hard_cor, start = "bari")$ns, 3)
})

test_that("FW with incorrect hard seeds", {
  expect_equal(graph_match_FW(A, B, seeds = hard_incor, start = "bari")$ns, 3)
})