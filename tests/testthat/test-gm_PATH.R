

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(10, .9, .5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
seeds<-1:4
tt <- graph_match_PATH(A, B, seeds)


test_that("matching correspondence between graph1 and graph2", {
  set.seed(123)
  cgnp_pair <- sample_correlated_gnp_pair(10, .9, .5)
  A <- cgnp_pair$graph1
  B <- cgnp_pair$graph2
  seeds <- 1:4
  tt <- graph_match_PATH(A, B, seeds)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
  
})


# startm <- matrix(rnorm(100), 10)
# test_that("add similarity scores", {
#   expect_snapshot_value(
#     graph_match_PATH(A, B, seeds, similarity = startm), 
#     "serialize"
#   )
# })

# # sample a pair of directed graphs 
# set.seed(123)
# cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr = .3, p = .5, directed = TRUE)
# A <- cgnp_pair$graph1
# B <- cgnp_pair$graph2
# seeds <- c(1,3,5)

# test_that("PATH for directed graphs", {
#   expect_snapshot_value(
#     graph_match_PATH(A, B, seeds),
#     "serialize"
#   )
# })


# set.seed(12)
# gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
# A <- lapply(gp_list, function(gp)gp[[1]])
# B <- lapply(gp_list, function(gp)gp[[2]])
# seeds <- 1:3

# test_that("PATH multi-layer", {
#   expect_snapshot_value(
#     graph_match_PATH(A, B, seeds = 1:3, epsilon = 5),
#     "serialize"
#   )
# })

