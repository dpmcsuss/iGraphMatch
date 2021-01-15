context("lap")

set.seed(12345)
cost <- Matrix::rsparsematrix(10, 10, .5)

a1 <- do_lap(cost, "lapjv")
a2 <- do_lap(cost, "lapmod")
a3 <- do_lap(cost, "clue")

sol <- c(1, 7, 9, 8, 3, 5, 6, 4, 10, 2)
test_that(
  "lapjv works",
  {
    expect_equivalent(a1, sol)
  }
)

test_that(
  "lapmod works",
  {
    expect_equivalent(a2, sol)
  }
)

test_that(
  "clue works",
  {
    expect_equivalent(a3, sol)
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
A <- cgnp_pair$graph1
B <- igraph::induced_subgraph(cgnp_pair$graph2, 1:10)

test_that(
  "error on other",
  {
    expect_error({
      graph_match_FW(A, B, start = "bari", lap_method = "other")
    },
    "Unrecognized LAP method: other.*")
  }
)