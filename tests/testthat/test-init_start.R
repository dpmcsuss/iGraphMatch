
set.seed(123)

# not testing convex b/c of how janky it is

set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 5, corr =  0.5, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
hs <- 1:5 <= 2
ss <- t(matrix(c(3, 4)))
suppressWarnings(
  res <- init_start(start = "convex", nns = 3, ns=2,
    soft_seeds = ss, A = g1[], B = g2[], seeds = hs)
)


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
  expect_snapshot_output(
    init_start(start = "rds", nns = 2)
  )
})
set.seed(123)
test_that("doubly stochastic matrix start w/o soft seeds",
  {
    expect_snapshot_output(
      init_start(start = "rds_perm_bari", nns = 2)
    )
  }
)


# initialize start matrix with soft seeds

ss<-t(matrix(c(3,4)))
set.seed(123)
test_that("bari start w. soft seeds", {
  start <- init_start(start = "bari", nns = 3,ns=2,soft_seeds=ss)
  expect_snapshot_output(start)
  expect_s4_class(start, "splrMatrix")
})

set.seed(123)
test_that("random doubly stochastic start w. soft seeds", {
  start <- init_start(start = "rds", nns = 3,ns=2,soft_seeds=ss)
  expect_snapshot_output(start)
})

set.seed(123)
test_that("doubly stochastic matrix start w. soft seeds",
  {
    expect_snapshot_output(
        init_start(start = "rds_perm_bari", nns = 3,ns=2,soft_seeds=ss)
    )
  }
)

# not testing convex b/c of how janky it is

set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 5, corr =  0.5, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
hs <- 1:5 <= 2
ss <- t(matrix(c(3, 4)))
suppressWarnings(
  res <- init_start(start = "convex", nns = 3, ns=2,
    soft_seeds = ss, A = g1[], B = g2[], seeds = hs)
)

expected <- structure(
  c(0.786, 0.053, 0.081,
    0.161, 0.786, 0.053,
    0.053, 0.161, 0.786),
  .Dim = c(3L, 3L), .Dimnames = list(NULL, NULL))
test_that("convex start w. soft seeds", {
  expect_snapshot_output(res)
})




test_that(
  "Error on overspecified soft seeds",
  {
    expect_error(
      {
        init_start(matrix(1, 4, 4), 4, soft_seeds = c(1,3))
      },
      "You are trying to use soft seeds but .*"
    )
  }
)


f <- function(nns,ns, soft_seeds) {
  matrix(0, nns, nns)
}
test_that(
  "Function as start",{
    expect_snapshot_output(
      init_start(f, 10)
    )
  }
)

f <- function(){}
test_that(
  "Function as start wrong args",
  {
    expect_error(
      init_start(f, 10),
      ".*functions passed to init_start must have at least the arguments nns, ns, and softs_seeds"
    )
  }
)


f <- function(nns, ns, soft_seeds) {
  n <- nns - nrow(check_seeds(soft_seeds, nns + ns)$seeds)
  matrix(runif((n - 1) ^ 2), n - 1)
}
test_that(
  "Function as start wrong size",
  {
    expect_error(
      init_start(f, 10),
      ".*must return a square matrix-like object with dimension"
    )
  }
)



test_that(
  "Invalid string",
  {
    expect_error(
      init_start("string", 10),
      "start must be either a matrix, function, or one of.*"
    )
  }
)


test_that(
  "RDS from sim start",
  {
    sim <- matrix(runif(9), 3)
    m <- init_start(start = "rds_from_sim", nns = 3, sim = sim)
    expect_snapshot_output(print(round(m, 4)))


    sim <- Matrix::rsparsematrix(10, 10, .4, rand.x = function(n) rep(1,n))
    m <- init_start(start = "rds_from_sim", nns = 3, sim = sim)
    expect_snapshot_output(print(round(m, 4)))

    expect_error(
      init_start(start = "rds_from_sim", nns = 3, sim = "asdf"),
      "Error: sim must be a matrix-like object.*"
    )

  }
)

test_that(
  "soft seeds with non-initial seeds",
  {
    seeds <- c(1, 4)
    soft_seeds <- c(2, 3)

    s <- init_start("bari", nns = 5, seeds = seeds, soft_seeds = soft_seeds)

    expect_equal(diag(s@x), c(1, 1, 0, 0, 0))
  }
)