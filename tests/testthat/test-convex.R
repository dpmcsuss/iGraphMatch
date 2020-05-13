context("Convex")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n=10, corr=0.3, p=0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
seeds <- 1:10 <= 3

ex_df <- data.frame(corr_A = c(1:10),
          corr_B = c(4, 2, 9, 8, 5, 7, 10, 6, 3, 1))

actual <- graph_match_convex(A, B)

test_that("matching correspondence between graph1 and graph2",
  {
    expect_equal(actual$corr, ex_df)
  })
test_that("permutation matrix", {
  expect_equal(actual$P,
    get_perm(10, 10, ex_df))
})
test_that("number of seeds", {
  expect_equal(actual$ns,0)
})




