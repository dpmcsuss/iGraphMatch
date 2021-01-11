context("hard-seeding Frank Wolfe with incorrect seeds")

# sample pair of graphs w. 10 vertices
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2

# sample hard seeds
hard_cor <- c(1, 5, 3) # correct
hard_incor <- data.frame(A = c(1, 2, 3), B = c(1, 2, 5)) # contain incorrect

test_that("FW with no seed", {
  expect_equal(nrow(graph_match_FW(A, B, seeds = NULL, start = "bari")$seeds), 0)
})

test_that("FW with correct hard seeds", {
  expect_equal(graph_match_FW(A, B, seeds = hard_cor, start = "bari")$seeds, data.frame(A = c(1, 5, 3), B = c(1, 5, 3)))
})

test_that("FW with incorrect hard seeds", {
  expect_equal(graph_match_FW(A, B, seeds = hard_incor, start = "bari")$seeds, data.frame(A = c(1, 2, 3), B = c(1, 2, 5)))
})

set.seed(12)
g2 <- sample_correlated_gnp_pair(n = 10, corr = 0.7, p = 0.3)
A2 <- g2$graph1
B2 <- g2$graph2

A_l <- list(A, A2)
B_l <- list(B, B2)

test_that("FW multi-layer", {
  expect_equal(graph_match_FW(A_l, B_l, seeds = 1:3, start = "bari")$seeds, data.frame(A = 1:3, B = 1:3))
})


# test get_s_to_ns, input non-list adj matrices 
seeds <- check_seeds(seeds = 1:3, nv = 10)
nonseeds <- seeds$nonseeds
seeds <- seeds$seeds
test_that("get_s_to_ns non-list adj matrices", {
  expect_equal(nrow(get_s_to_ns(A, B, seeds, nonseeds, perm = sample(10-nrow(seeds)))), 7)
})


