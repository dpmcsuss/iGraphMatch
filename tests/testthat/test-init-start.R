context("initialization of the start matrix")


## initialize start matrix without soft seeds
set.seed(123)
test_that("bari start w/o soft seeds", {
  expect_equal(init_start(start = "bari", nns = 2),
               matrix(c(0.5,0.5,0.5,0.5),nrow=2))
})
set.seed(123)
test_that("random doubly stochastic start w/o soft seeds", {
  expect_equal(round(init_start(start = "rds", nns = 2),2),
               matrix(c(0.47,0.53,0.53,0.47),nrow=2))
})
set.seed(123)
test_that("doubly stochastic matrix start w/o soft seeds", {
  expect_equivalent(round(as.matrix(init_start(start = "rds_perm_bari", nns = 2)),2),
               matrix(c(0.36,0.64,0.64,0.36),nrow=2))
})


# initialize start matrix with soft seeds












