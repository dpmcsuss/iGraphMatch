context("best matches ranking")

# sample pair of graphs w. 10 vertices
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2
seeds <- 1:10 <= 3
nonseeds <- !seeds
match_corr <- graph_match_FW(A, B, seeds, start = "bari")$corr

test_that("best matches for all matches using all measures", {
  expect_equal(nrow(best_matches(A, B, measure = "row_cor", 
                                         num = 3, x = nonseeds, match_corr$corr_B)), 3)
  expect_equal(nrow(best_matches(A, B, measure = "row_diff", 
                                         num = 3, x = nonseeds, match_corr$corr_B)), 3)
  expect_equal(nrow(best_matches(A, B, measure = "row_perm_stat", 
                                         num = 3, x = nonseeds, match_corr$corr_B)), 3)
})
