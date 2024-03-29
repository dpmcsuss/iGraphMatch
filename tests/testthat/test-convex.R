

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n=10, corr=0.8, p=0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2

# ex_df <- data.frame(corr_A = c(1:10),
          # corr_B = c(4, 2, 9, 8, 5, 7, 10, 6, 3, 1))

test_that("warning from gm convex",
  {
    expect_warning(
      actual <<- gm(A, B, method = "convex"),
      "Frank-Wolfe iterations reach the maximum iteration, convergence may not occur.*"
    )
  }
)

# test_that("correct matching result",
#   {

#     expect_snapshot_value(actual@corr, style = "serialize")
#   })

test_that("matching correspondence between graph1 and graph2",
  {
    expect_equal(dim(actual), c(igraph::vcount(A), igraph::vcount(B)))
  }
)
test_that("doubly stochastic", {
  expect_lt(sum(abs(rowSums(actual$soft) - 1)), 10e-6)
  expect_lt(sum(abs(colSums(actual$soft) - 1)), 10e-6)
})
test_that("number of seeds", {
  expect_equal(sum(actual$seeds), 0)
})


# test output error when given start = "convex"
test_that("doubly stochastic", {
  expect_error(
    gm(A, B, method = "convex",start = "convex"),
    "Cannot start convex with convex. Try \"bari\" or another option."
  )
})
