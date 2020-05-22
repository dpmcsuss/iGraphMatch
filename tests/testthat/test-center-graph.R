context("center graph")

set.seed(123)
A <- sample_correlated_gnp_pair(n = 3, corr = .5, p = .5)$graph1
test_that("naive scheme", {
  expect_equivalent(as.matrix(center_graph(A, scheme = "naive")),
               matrix(c(0,0,1,0,0,0,1,0,0),nrow=3))
})
test_that("center scheme", {
  expect_equivalent(as.matrix(center_graph(A, scheme = "center")),
               matrix(c(-1,-1,1,-1,-1,-1,1,-1,-1),nrow=3))
})
test_that("scheme equal to 1", {
  expect_equivalent(as.matrix(center_graph(A, scheme = 1)),
               matrix(c(0,0,1,0,0,0,0,0,0),nrow=3))
})
