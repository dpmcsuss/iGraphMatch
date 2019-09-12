context("measure functions")

cgnp_pair <- sample_correlated_gnp_pair(n = 50, corr =  0.3,p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
match <- graph_match_FW(g1, g2)
g2m <- g2[match$corr$corr_B, match$corr$corr_B]
g1 <- g1[]

test_that("measure functions", {
  expect_vector(row_cor(g1, g2m), size = 50)
  expect_vector(row_diff(g1, g2m), size = 50)
  expect_vector(row_perm_stat(g1, g2m), size = 50)
})


