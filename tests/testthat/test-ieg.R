context("ieg")

#Sample graphs from edge probability matrix and correlation matrix
set.seed(123)
n <- 10
p_mat <- matrix(runif(n^2),n)
c_mat <- matrix(runif(n^2),n)
g<-sample_correlated_ieg_pair(n, p_mat, c_mat, directed = TRUE,loops = TRUE, permutation = 1:n)

test_that("number of vertices", {
  expect_equal(length(V(g$graph1)),10) 
  expect_equal(length(V(g$graph2)),10)
})
test_that("number of edges", {
  expect_equal(length(E(g$graph1)),56)
  expect_equal(length(E(g$graph2)),56)
})
test_that("degree of vertex in each graph", {
  expect_equal(degree(g$graph1),c(10,10,13,13,11,9,13,13,12,8))
  expect_equal(degree(g$graph2),c(11,10,10,11,9,13,13,13,11,11)) 
})