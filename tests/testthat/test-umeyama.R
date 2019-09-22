context("Umeyama")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(10, .9, .5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
startm <- matrix(0, 10, 10)
diag(startm)[1:4] <- 1


test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_Umeyama(A, B,startm, alpha = .3)$corr,
               data.frame(corr_A = c(1:10), corr_B = c(1,2,3,4,5,10,7,8,9,6)))
})
test_that("number of seeds", {
  expect_equal(graph_match_Umeyama(A, B,startm, alpha = .3)$ns,0)
})


