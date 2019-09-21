context("IsoRank")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
A <- g$graph1
B <- g$graph2
startm <- matrix(0, 10, 10)
diag(startm)[1:4] <- 1


test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_IsoRank(A, B,startm, alpha = .3, method = "greedy")$order,
               data.frame(corr_A = c(3,4,1,2,9,5,8,7,10,6), corr_B = c(3,4,1,2,5,9,8,10,7,6)))
})
test_that("order of nodes getting matched", {
  expect_equal(graph_match_IsoRank(A, B,startm, alpha = .3, method = "greedy")$order,
               c(3,4,1,2,6,10,8,7,5,9))
})
test_that("test number of seeds", {
  expect_equal(graph_match_IsoRank(A, B,startm, alpha = .3, method = "greedy")$ns,0)
})