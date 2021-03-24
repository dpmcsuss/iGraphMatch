

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(10, .9, .5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
startm <- matrix(0, 10, 10)
diag(startm)[1:4] <- 1
seeds<-1:4




test_that("matching correspondence between graph1 and graph2", {
  tt <- graph_match_Umeyama(A, B, seeds, startm)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")

})
# test_that("number of seeds", {
#   expect_equal(graph_match_Umeyama(A, B, seeds, startm)$seeds, 
#                data.frame(A = 1:4, B = 1:4))
# })

# sample a pair of directed graphs 
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr = .9, p = .5, directed = TRUE)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2

test_that("matching correspondence between graph1 and graph2 for directed graphs", {
  tt <- graph_match_Umeyama(A, B, seeds, startm)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
})

# test_that("number of seeds for directed graphs", {
#   expect_equal(nrow(graph_match_Umeyama(A, B, similarity = startm)$seeds), 0)
# })


set.seed(12)
gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
B <- lapply(gp_list, function(gp)gp[[2]])
seeds <- 1:3

test_that("Umeyama multi-layer", {
  tt <- graph_match_Umeyama(A, B, seeds)
  expect_snapshot_output(print(tt))
  expect_snapshot_value(tt, "serialize")
})


