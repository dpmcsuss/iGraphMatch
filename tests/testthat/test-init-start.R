context("initialization of the start matrix")


## initialize start matrix without soft seeds
set.seed(123)
test_that("bari start w/o soft seeds", {
  expected <- splr(Matrix(0,2,2),
    a = Matrix(c(1,1),2), b = Matrix(c(.5,.5),2))
  expect_equal(init_start(start = "bari", nns = 2),
               expected)
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

ss<-t(matrix(c(3,4)))
set.seed(123)
test_that("bari start w. soft seeds", {
  expect_equivalent(init_start(start = "bari", nns = 3,ns=2,soft_seeds=ss),
               matrix(c(0,0.5,0.5,1,0,0,0,0.5,0.5)))
})

set.seed(123)
test_that("random doubly stochastic start w. soft seeds", {
  expect_equivalent(round(init_start(start = "rds", nns = 3,ns=2,soft_seeds=ss),2),
               matrix(c(0,0.47,0.53,1,0,0,0,0.53,0.47)))
})

set.seed(123)
test_that("doubly stochastic matrix start w. soft seeds", {
  expect_equivalent(round(as.matrix(init_start(start = "rds_perm_bari", nns = 3,ns=2,soft_seeds=ss)),2),
                    matrix(c(0,0.36,0.64,1,0,0,0,0.64,0.36)))
})

set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 5, corr =  0.5, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
hs <- 1:5 <= 2
ss <- t(matrix(c(3, 4)))
res <- init_start(start = "convex", nns = 3, ns=2, 
      soft_seeds = ss, A = g1, B = g2, seeds = hs)
test_that("doubly stochastic matrix start w. soft seeds", {
  expect_equivalent(
    as.matrix(res),
    t(matrix(c(
          0, 1,   0,
        0.7, 0, 0.3,
        0.3, 0, 0.7), 3)))
})
