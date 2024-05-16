

start <- init_start("bari", 2, 2, soft_seeds = data.frame(A = 2, B = 2))
test_that("bari start w. soft seeds", {
  expect_equal(as.matrix(splr_to_sparse(start)), diag(2))
})

test_that("doubly stochastic matrix start w. soft seeds", {
  expect_equal(
    as.matrix(round(init_start("rds_from_sim", 2, 2, sim = diag(1e6, 2, 2)))),
    diag(2)
  )
})

sim <- Matrix::rsparsematrix(10, 10, .4, rand.x = function(n) rep(1,n))
test_that("doubly stochastic matrix start w. soft seeds", {
  s <- rds_from_sim_start(10, sim = sim)
  expect_equal(nrow(s), 10)
  s <- rds_from_sim_start(10, sim = as.matrix(sim))
  expect_equal(nrow(s), 10)
  expect_snapshot_output({
    print(round(s, 4))
  })
  expect_warning(
    rds_from_sim_start(10, soft_seeds = data.frame(A = 4, B = 4), sim = sim),
    "Ignoring soft_seeds in rds_from_sim_start")
})

