context("measure functions")

set.seed(1234)
n <- 10
cgnp_pair <- sample_correlated_gnp_pair(n = n, corr =  0.3,p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
match <- graph_match_FW(g1, g2, start = "bari")
g2m <- g2[match$corr$corr_B, match$corr$corr_B]
g1 <- g1[]

test_that("measure functions", {
  expect_length(row_cor(g1, g2m), n)
  expect_length(row_diff(g1, g2m), n)
  expect_length(row_perm_stat(g1, g2m), n)
  expect_length(row_perm_stat(g1, g2m, exact = FALSE), n)
})

test_that("check sym case", {
  expect_length(row_diff_perm(g2, g2m, nmc = 100, sym = TRUE)$mean, n)
  expect_length(mean_row_diff(g2, g2m, sym = TRUE), n)
  expect_length(var_row_diff(g2, g2m, sym = TRUE), n)
})

