# largest common cc

set.seed(123)
n <- 10
cgnp_pair <- sample_correlated_gnp_pair(n = n,
  corr =  0.7, p =  0.2)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
lcc1 <- largest_common_cc(g1, g2, min_degree = 1)
lcc3 <- largest_common_cc(g1, g2, min_degree = 3)

test_that("largest cc w. min_degree", {
  expect_length(lcc1, 3)
  expect_equal(sum(lcc1$keep), 4)
  expect_length(lcc1$keep, n)

  expect_length(lcc3, 3)
  expect_equal(sum(lcc3$keep), 2)
  expect_length(lcc3$keep, n)
})


# largest cc

set.seed(123)
g <- igraph::sample_gnp(100, .01)
lcc <- largest_cc(g)

test_that("largest cc", {
  expect_length(lcc, 2)
  expect_equal(sum(lcc$keep), 46)
  expect_length(lcc$keep, 100)
})
