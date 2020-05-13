context("start initialization")

test_that("bari start w. soft seeds", {
  start <- bari_start(2, 2, 
    soft_seeds = data.frame(A = 4, B = 4))
  expect_equal(start, diag(2))
})

test_that("random doubly stochastic start w. soft seeds", {
    expect_equal(rds_sinkhorn_start(2, 2, soft_seeds = data.frame(A = 4, B = 4)),
               diag(2))
})

test_that("doubly stochastic matrix start w. soft seeds", {
  expect_equal(rds_sinkhorn_start(2, 2, soft_seeds = data.frame(A = 4, B = 4)),
               diag(2))
})

