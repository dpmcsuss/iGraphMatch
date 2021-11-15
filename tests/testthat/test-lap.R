

set.seed(12345)
cost <- Matrix::rsparsematrix(10, 10, .5)


sol <- c(1, 7, 9, 8, 3, 5, 6, 4, 10, 2)
test_that(
  "lapjv works",
  {
    a1 <- do_lap(cost, "lapjv") 
    expect_equal(a1, sol)
  }
)

test_that(
  "lapmod works",
  {
    a2 <- do_lap(cost, "lapmod")
    expect_equal(a2, sol)
  }
)

test_that(
  "clue works",
  {
    a3 <- do_lap(cost, "clue")
    expect_equal(a3, sol)
  }
)

test_that(
  "rect works",
  {
    a4 <- do_lap(cost, "rect")
    expect_equal(a4, sol)
  }
)

test_that(
  "error on other",
  {
    expect_error({
      do_lap(cost, "other")
    },
    "Unrecognized LAP method: other Please use one of: lapmod lapjv clue.*")
  }
)

set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n=100, corr=0.8, p=0.5)
B <- cgnp_pair$graph1
A <- igraph::induced_subgraph(cgnp_pair$graph2, 1:10)

test_that(
  "error on other",
  {
    expect_error({
      gm(A, B, start = "bari", method = "indefinite", lap_method = "other")
    },
    "Unrecognized LAP method: other.*")
  }
)


test_that(
  "harder lap",
  {
    expect_snapshot_output({
      hard <-
        matrix(sample(1000, 10000, replace = TRUE), 100) +
        sample(1000, 100, replace = TRUE)
      print(lapmod(hard))
    })
  }
)

