context("largest cc")

cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.7, p =  0.2)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2

test_that("largest cc w. min_degree", {
  expect_vector(largest_common_cc(g1, g2, min_degree = 1))
  expect_vector(largest_common_cc(g1, g2, min_degree = 3))
})