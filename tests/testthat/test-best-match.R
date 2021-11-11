

# sample pair of graphs w. 10 vertices
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.6, p = 0.5, directed = TRUE)
A <- g$graph1
B <- g$graph2
seeds <- 1:3
match <- gm(A, B, seeds, r = 3, method = "percolation")

test_that("best matches for all matches using all measures", {
  expect_equal(nrow(best_matches(A, B, match, measure = "row_cor",
                                         num = 4)), 4)
  expect_equal(nrow(best_matches(A, B, match, measure = "row_diff",
                                 num = 4)), 4)
  expect_equal(nrow(best_matches(A, B, match, measure = "row_perm_stat",
                                 num = 4)), 4)
})


test_that(
  "Error on bad measure",
  {
    expect_error(
      best_matches(A, B, match, measure = "hello", 4),
      "*measure must be one of*"
    )
    expect_error(
      best_matches(A, B, match, measure = function(a, b) c(a, b), 4),
      "*measure must be one of*"
    )
  }
)

rel_row_diff <- function(g1, g2){
  g1 <- g1[]
  g2 <- g2[]
  g1 <- as.matrix(g1)
  g2 <- as.matrix(g2)
  rowSums(abs(g1-g2)) / (rowSums(g1) + rowSums(g2))
}

test_that(
  "Function works for measure",
  {
    expect_equal(
      nrow(best_matches(A, B, match, measure = rel_row_diff, 4)),
      4
    )
  }
)

test_that(
  "true correspondence available",
  {
    expect_equal(ncol(best_matches(A, B, match, measure = "row_cor",
                                   num = 4, true_label = 1:10)), 4)
    expect_equal(ncol(best_matches(A, B, match, measure = "row_diff",
                                   num = 4, true_label = 1:10)), 4)
    expect_equal(ncol(best_matches(A, B, match, measure = "row_perm_stat",
                                   num = 4, true_label = 1:10)), 4)
  }
)
