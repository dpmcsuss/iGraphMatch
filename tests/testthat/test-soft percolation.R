context("Soft Percolation")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
seeds <- 1:5


test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_soft_percolation(A, B,seeds,r=2, max_iter = 2)$corr,
               data.frame(corr_A = c(1,2,3,4,5,9), corr_B = c(1,2,3,4,5,7)))
})
test_that("order of nodes getting matched", {
  expect_equal(graph_match_soft_percolation(A, B,seeds,r=2, max_iter = 2)$order,
               c(1,2,3,4,5,6))
})
test_that("test number of seeds", {
  expect_equal(graph_match_soft_percolation(A, B,seeds,r=2, max_iter = 2)$ns,5)
})