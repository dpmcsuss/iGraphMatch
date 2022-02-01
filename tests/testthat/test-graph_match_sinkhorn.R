
cgnp_pair <- sample_correlated_gnp_pair(n = 20, corr =  0.9, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
# match G_1 & G_2 with no seeds

test_that("gm_sinkhorn", {
  expect_snapshot_output(
    print(
      suppressWarnings(
        gm(g1, g2, method = "sinkhorn", max_iter = 10)
      )
    )
  )
  seeds <- 1:20 <= 3
  expect_snapshot_output(
    print(
      suppressWarnings(
        gm(g1, g2, seeds, method = "sinkhorn", max_iter = 10)
      )
    )
  )
})