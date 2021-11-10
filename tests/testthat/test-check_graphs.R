

# sample pair of graphs w. 10 vertices
cgnp_pair <- sample_correlated_gnp_pair(n=10, corr=0.8, p=0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2

g <- check_graph(A, B, as_igraph = TRUE)


test_that(
  "Correct classes",
  {
    expect_s3_class(g$g1, "igraph")
    expect_s3_class(g$g2, "igraph")
  }
)
test_that(
  "Error on not igraph",
  {
    expect_error(
      {
        check_graph(A[], B, as_igraph = TRUE)
      },
      "Check graph only supports as_igraph = TRUE if both A and B are igraph objects."
    )
  }
)


check_graph(igraph::induced_subgraph(A, 1:4), B, as_igraph = TRUE)
check_graph(A, igraph::induced_subgraph(B, 1:4), as_igraph = TRUE)

Al <- list(matrix(runif(4^2), 4), matrix(runif(5^2),5))
Bl <- list(matrix(runif(5^2), 5), matrix(runif(5^2),5))

test_that(
  "Error on not equal order",
  {
    expect_error(
      {
        check_graph(Al, Bl)
      },
      "A contains graphs of different orders. For multiple graph matching, all graphs must have the same number of vertices."
    )
    expect_error(
      {
        check_graph(Bl, Al)
      },
      "B contains graphs of different orders. For multiple graph matching, all graphs must have the same number of vertices."
    )
  }
)

Al <- list(matrix(runif(4^2), 4), matrix(runif(4^2), 4))
check_graph(Al[[1]], Bl[[1]])
check_graph(Al[[1]], Bl[[1]], as_list = FALSE)
check_graph(Al, Bl)
check_graph(Bl, Al)


test_that(
  "Error not as_list for multilayer",
  {
    expect_error(
      {
        check_graph(Al, Bl, as_list = FALSE)
      },
      "A is multi-layer and must be converted to single layer.*"
    )
    expect_error(
      {
        check_graph(Al[[1]], Bl, as_list = FALSE)
      },
      "B is multi-layer and must be converted to single layer.*"
    )
  }
)






Al <- list(matrix(runif(4^2), 8), matrix(runif(4^2), 8))



test_that(
  "Error not square",
  {
    expect_error(
      {
        check_graph(Al, Bl)
      },
      "A is not square. This method only supports.*"
    )
    expect_error(
      {
        check_graph(Bl, Al)
      },
      "B is not square. This method only supports.*"
    )
    expect_error(
      {
        check_graph(Al, Bl, square = FALSE)
      },
      "square = FALSE is not yet supported for check_graph"
    )
  }
)
