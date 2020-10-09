context("ExpandWhenStuck")

set.seed(12)
G <- sample_correlated_gnp_pair(n = 10, corr = .5, p = .5)
A <- G$graph1
B <- G$graph2
seeds <- 1:4

test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds, r = 2)$corr,
               data.frame(corr_A = c(1:7, 9:10), corr_B = c(1:4, 9, 6, 7, 5, 10),row.names=as.character(1:9)))
})

test_that("test number of seeds", {
  expect_equal(graph_match_ExpandWhenStuck(A, B,seeds, r = 2)$ns, 4)
})

# with similarity score
sim <- matrix(rnorm(100), 10)
test_that("perco w. similarity score", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds, similarity = sim)$ns, 4)
})

# directed graphs
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5, directed = TRUE)
A <- g$graph1
B <- g$graph2
test_that("perco w. similarity score", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds, similarity = sim)$ns, 4)
})

# multi-layer
set.seed(1)
g2 <- sample_correlated_gnp_pair(n = 10, corr = 0.7, p = 0.3)
A2 <- g2$graph1
B2 <- g2$graph2

A_l <- list(A, A2)
B_l <- list(B, B2)

test_that("ExpandWhenStuck multi-layer", {
  expect_equal(graph_match_ExpandWhenStuck(A_l, B_l, seeds = 1:3)$ns, 3)
})


# multiple candidate matches with same score
set.seed(123)
g <- sample_correlated_gnp_pair(n = 50, corr = 0.5, p = 0.8, directed = TRUE)
A <- g$graph1
B <- g$graph2
seeds <- 1:10
test_that("perco w. similarity score", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds)$ns, 10)
})



