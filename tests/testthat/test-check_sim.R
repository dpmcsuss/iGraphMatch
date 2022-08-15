

set.seed(123)
n <- 20
cgnp_pair <- sample_correlated_gnp_pair(n=n, corr=0.8, p=0.5)
A <- cgnp_pair$graph1
B <- igraph::induced_subgraph(cgnp_pair$graph2, 1:10)

sim_rect <- matrix(runif(10 * n), n)
sim_sq <- pad(sim_rect, 0, n - 10)
gm(A, B, start = "bari", sim = sim_rect, method = "indefinite")
gm(A, B, start = "bari", sim = sim_sq, method = "indefinite")
gm(A, B, start = "bari", seeds = 1:3, sim = sim_sq, method = "indefinite")
gm(A, B, start = "bari", sim = NULL, method = "indefinite")


test_that(
  "error on wrong rectangular dim",
  {
    expect_error({
      sim_bad <- matrix(runif(10 * n - 10), n - 10)
      gm(A, B, start = "bari", sim = sim_bad, method = "indefinite")
    },
    "Non square similarity matrices.*")
  }
)


test_that(
  "error on wrong square dim",
  {
    expect_error({
      sim_bad <- matrix(runif(10 * 10), 10)
      gm(A, B, start = "bari", sim = sim_bad, method = "indefinite")
    },
    "Square similarity matrices.*")
  }
)
