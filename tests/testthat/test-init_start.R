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
  expect_equivalent(
    as.matrix(round(init_start(start = "rds", nns = 2),2)),
    matrix(c(0.47,0.53,0.53,0.47),nrow=2))
})
set.seed(123)
test_that("doubly stochastic matrix start w/o soft seeds", 
  {
    expect_equivalent(
      round(as.matrix(init_start(start = "rds_perm_bari", nns = 2)),2),
      structure(
        c(0.64, 0.36, 0.36, 0.64),
        .Dim = c(2L, 2L), .Dimnames = list(NULL, NULL))
    )
  }
)


# initialize start matrix with soft seeds

ss<-t(matrix(c(3,4)))
set.seed(123)
test_that("bari start w. soft seeds", {
  start <- init_start(start = "bari", nns = 3,ns=2,soft_seeds=ss)
  expect_equivalent(as.matrix(start),
               matrix(c(0,0.5,0.5,1,0,0,0,0.5,0.5)))
  expect_s4_class(start, "splrMatrix")
})

set.seed(123)
test_that("random doubly stochastic start w. soft seeds", {
  start <- init_start(start = "rds", nns = 3,ns=2,soft_seeds=ss)
  expect_equivalent(round(as.matrix(start), 2),
               matrix(c(0,0.47,0.53,1,0,0,0,0.53,0.47)))
})

set.seed(123)
test_that("doubly stochastic matrix start w. soft seeds", 
  {
    expect_equivalent(
      round(as.matrix(
        init_start(start = "rds_perm_bari", nns = 3,ns=2,soft_seeds=ss)
      ),2),
      structure(
        c(0, 0.64, 0.36, 1, 0, 0, 0, 0.36, 0.64),
        .Dim = c(3L, 3L), .Dimnames = list(NULL, NULL)
      )
    )
  }
)


set.seed(123)

# not testing convex b/c of how janky it is

set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 5, corr =  0.5, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
hs <- 1:5 <= 2
ss <- t(matrix(c(3, 4)))
res <- init_start(start = "convex", nns = 3, ns=2, 
      soft_seeds = ss, A = g1[], B = g2[], seeds = hs)

expected <- structure(
  c(
    0.668324233550846, 0.100637887848042, 0.231037878601112, 
    0.135904634471242, 0.806794628456265, 0.0573007370724925,
    0.195771131977911, 0.0925674836956933, 0.711661384326395
  ),
  .Dim = c(3L, 3L), .Dimnames = list(NULL, NULL))
test_that("convex start w. soft seeds", {
  expect_equivalent(
    as.matrix(res),
    expected
  )
})




test_that(
  "Error on overspecified soft seeds",
  expect_error(
    {
      init_start(matrix(1, 4, 4), 4, soft_seeds = c(1,3))
    },
    "You are trying to use soft seeds but .*"
  )
)


f <- function(nns,ns, soft_seeds) {
  matrix(0, nns, nns)
}
test_that(
  "Function as start",
  expect_equivalent(
    {
      dim(init_start(f, 10))
    },
    c(10, 10)
  )
)

f <- function(){}
test_that(
  "Function as start wrong args",
  expect_error(
    {
      init_start(f, 10)
    },
    ".*functions passed to init_start must have at least the arguments nns, ns, and softs_seeds"
  )
)


f <- function(nns, ns, soft_seeds) {
  n <- nns - nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  matrix(runif((n - 1) ^ 2), n - 1)
}
test_that(
  "Function as start wrong size",
  expect_error(
    {
      init_start(f, 10)
    },
    ".*must return a square matrix-like object with dimension"
  )
)



test_that(
  "Invalid string",
  expect_error(
    {
      init_start("string", 10)
    },
    "Start must be either a matrix, function, or one of.*"
  )
)
