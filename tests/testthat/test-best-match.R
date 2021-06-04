

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
