

# sample pair of graphs w. 10 vertices
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2
seeds <- c(1, 5, 3)

test_that("perco of same sizes", {
  tt <- gm(A, B, seeds, method = "percolation", ExpandWhenStuck = FALSE)
  expect_snapshot_output(tt)
  expect_snapshot_value(tt, "serialize")
})

# with similarity score
sim <- matrix(rnorm(100), 10)
test_that("perco w. similarity score", {
  tt <- gm(A, B, seeds, similarity = sim, method = "percolation", ExpandWhenStuck = FALSE)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
})

test_that("percolation without seeds", {
  tt <- gm(A, B, seeds = NULL, similarity = sim, method = "percolation", ExpandWhenStuck = FALSE)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
})


test_that("perco w. directed graphs", {
  # directed graphs
  set.seed(123)
  g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5, directed = TRUE)
  A <- g$graph1
  B <- g$graph2
  tt <- gm(A, B, seeds, similarity = sim, method = "percolation", ExpandWhenStuck = FALSE)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
})



test_that("percolation multi-layer", {
  set.seed(12)
  gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
  A <- lapply(gp_list, function(gp)gp[[1]])
  B <- lapply(gp_list, function(gp)gp[[2]])
  seeds <- 1:3
  tt <- gm(A, B, seeds, method = "percolation", ExpandWhenStuck = FALSE)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
})

