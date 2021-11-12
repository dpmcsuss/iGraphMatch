

set.seed(123)
n <- 20
cgnp_pair <- sample_correlated_gnp_pair(n=n, corr=0.8, p=0.5)
B <- cgnp_pair$graph1
A <- igraph::induced_subgraph(cgnp_pair$graph2, 1:(n / 2))

sim_rect <- matrix(runif((n / 2) * n), n / 2)
sim_sq <- pad(sim_rect, n - (n / 2), 0)
gm(A, B, start = "bari", sim = sim_rect, method = "indefinite")
gm(A, B, start = "bari", sim = sim_sq, method = "indefinite")
gm(A, B, start = "bari", seeds = 1:3, sim = sim_sq, method = "indefinite")
gm(A, B, start = "bari", sim = NULL, method = "indefinite")


test_that(
  "error on wrong rectangular dim",
  {
    expect_error({
      sim_bad <- matrix(runif((n / 2) * n - (n / 2)), n - (n / 2))
      gm(A, B, start = "bari", sim = sim_bad, method = "indefinite")
    },
    ".*")
    # "Non square similarity matrices.*")
  }
)


test_that(
  "error on wrong square dim",
  {
    expect_error({
      sim_bad <- matrix(runif((n / 2) * (n / 2)), (n / 2))
      gm(A, B, start = "bari", sim = sim_bad, method = "indefinite")
    },
    ".*")
    # "Square similarity matrices.*")
  }
)
