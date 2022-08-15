

# sample pair of graphs w. 10 vertices
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2


test_that(
  "Graph match sinkhorn seems to work",
  {
    expect_snapshot_output(
      gm(A, B, method = "sinkhorn")$corr
    )
  }
)