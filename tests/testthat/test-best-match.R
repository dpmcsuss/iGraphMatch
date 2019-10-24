context("best matches ranking")

# sample pair of graphs w. 10 vertices
g <- sample_correlated_gnp_pair(n = 10, corr = 0.4, p = 0.5)
A <- g$graph1
B <- g$graph2
seeds <- 1:10 <= 3
nonseeds <- !seeds
match_corr <- graph_match_FW(A, B, seeds, start = "bari")$corr

test_that("best matches for all matches using all measures", {
  expect_equal(nrow(best_matches(A, B, measure = "row_cor", 
                                         num = 3, x = nonseeds, match_corr)), 3)
  expect_equal(nrow(best_matches(A, B, measure = "row_diff", 
                                         num = 3, x = nonseeds, match_corr)), 3)
  expect_equal(nrow(best_matches(A, B, measure = "row_perm_stat", 
                                         num = 3, x = nonseeds, match_corr)), 3)
})

cor_true <- data.frame(A_best = c(5,6,7), B_best = c(9,6,7))
diff_true <- data.frame(A_best = c(5,4,6), B_best = c(9,10,6))
perm_true <- cor_true
test_that("best matches for partial match result",{
  expect_equal(nrow(best_matches(A, B, measure = "row_cor", num = 3, x = nonseeds, 
                            match_corr[1:8,])), 3)
  expect_equal(nrow(best_matches(A, B, measure = "row_diff", num = 3, x = nonseeds, 
                            match_corr[1:8,])), 3)
  expect_equal(nrow(best_matches(A, B, measure = "row_perm_stat", num = 3, x = nonseeds, 
                            match_corr[1:8,])), 3)
})
