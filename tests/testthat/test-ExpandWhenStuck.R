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
  expect_equal(graph_match_ExpandWhenStuck(A, B,seeds, r = 2)$seeds,
               data.frame(A = 1:4, B = 1:4))
})

# with similarity score
sim <- matrix(rnorm(100), 10)
test_that("exp w. similarity score", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds, similarity = sim)$seeds,
               data.frame(A = 1:4, B = 1:4))
})

test_that("exp w. similarity score & no seeds", {
  expect_equal(nrow(graph_match_ExpandWhenStuck(A, B, seeds = NULL, similarity = sim)$seeds),
               0)
})

# directed graphs
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5, directed = TRUE)
A <- g$graph1
B <- g$graph2
test_that("perco w. similarity score", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds, similarity = sim)$seeds,
               data.frame(A = 1:4, B = 1:4))
})


# multiple candidate matches with same score
set.seed(12)
gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
B <- lapply(gp_list, function(gp)gp[[2]])
seeds <- 1:3
test_that("perco w. similarity score", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds = seeds)$seeds,
               data.frame(A = 1:3, B = 1:3))
})



