

# sample pair of graphs w. 10 vertices
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5)
A <- g$graph1
B <- g$graph2


# customized graph matching algorithm
graph_match_rand <- function(A, B, seeds = NULL, similarity = NULL, rand_seed){
  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  nv <- max(totv1, totv2)

  corr_A <- 1:nv
  set.seed(rand_seed)
  corr_B <- c(1:nv)[sample(nv)]
  corr <- data.frame(corr_A, corr_B)

  cl <- match.call()
  graphMatch(
    corr = corr,
    nnodes = c(totv1, totv2),
    call = cl,
    detail = list(
      rand_seed = rand_seed
    )
  )
}

test_that("customized gm method", {
  m_self <- gm(A, B, method = graph_match_rand, rand_seed = 123)
  expect_snapshot_output(m_self)
  expect_snapshot_value(m_self, "serialize")
})


test_that(
  "Error on invalid customized method",
  {
    expect_error(
      gm(A, B, method = "hello"),
      "*Method must be one of*"
    )
  }
)

test_that(
  "Error on invalid customized method output",
  {
    expect_error(
      gm(A, B, method = function(A, B, seeds, similarity) list(A, B)),
      "*Customized graph matching method function*"
    )
  }
)
