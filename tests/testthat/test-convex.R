context("Convex")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n=10, corr=0.8, p=0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2

# ex_df <- data.frame(corr_A = c(1:10),
          # corr_B = c(4, 2, 9, 8, 5, 7, 10, 6, 3, 1))

actual <- graph_match_convex(A, B)

test_that("matching correspondence between graph1 and graph2",
  {
    expect_equal(length(actual), 6)
  })
test_that("doubly stochastic", {
  expect_lt(sum(abs(rowSums(actual$D) - 1)), 10e-6)
  expect_lt(sum(abs(colSums(actual$D) - 1)), 10e-6)
})
test_that("number of seeds", {
  expect_equal(actual$ns,0)
})