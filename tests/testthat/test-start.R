context("start initialization")

test_that("bari start w. soft seeds", {
  start <- bari_start(2, 2, 
    soft_seeds = data.frame(A = 4, B = 4))
  expect_equivalent(as.matrix(start), diag(2))
})

test_that("random doubly stochastic start w. soft seeds", {
    expect_equivalent(
      as.matrix(rds_sinkhorn_start(2, 2, soft_seeds = data.frame(A = 4, B = 4))),
               diag(2))
})

test_that("doubly stochastic matrix start w. soft seeds", {
  expect_equivalent(
    as.matrix(rds_sinkhorn_start(2, 2, soft_seeds = data.frame(A = 4, B = 4))),
               diag(2))
})

sim <- Matrix::rsparsematrix(10, 10, .4, rand.x = function(n) rep(1,n))
test_that("doubly stochastic matrix start w. soft seeds", {
  expect_equal(
    nrow(rds_from_sim_start(10, sim = sim)),
    10)
  expect_equal(
    nrow(rds_from_sim_start(10, sim = as.matrix(sim))),
    10)
  expect_warning(
    nrow(rds_from_sim_start(10, soft_seeds = data.frame(A = 4, B = 4), sim = sim)),
    "Ignoring soft_seeds in rds_from_sim_start")
})

