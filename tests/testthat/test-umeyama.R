context("Umeyama")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(10, .9, .5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
startm <- matrix(0, 10, 10)
diag(startm)[1:4] <- 1
seeds<-1:4

test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_Umeyama(A, B, startm, seeds)$corr,
               data.frame(corr_A = c(1:10), corr_B = c(1,2,3,4,5,10,7,8,9,6)))
})
test_that("number of seeds", {
  expect_equal(graph_match_Umeyama(A, B, startm, seeds)$ns, 4)
})

# sample a pair of directed graphs 
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr = .9, p = .5, directed = TRUE)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2

test_that("matching correspondence between graph1 and graph2 for directed graphs", {
  expect_equal(graph_match_Umeyama(A, B, startm, seeds)$corr,
               data.frame(corr_A = c(1:10), corr_B = c(1,2,3,4,8,9,7,10,6,5)))
})
test_that("number of seeds for directed graphs", {
  expect_equal(graph_match_Umeyama(A, B, startm)$ns, 0)
})


set.seed(12)
g2 <- sample_correlated_gnp_pair(n = 10, corr = 0.7, p = 0.3)
A2 <- g2$graph1
B2 <- g2$graph2

A_l <- list(A, A2)
B_l <- list(B, B2)

test_that("Umeyama multi-layer", {
  expect_equal(graph_match_Umeyama(A_l, B_l, seeds = 1:3)$ns, 3)
})


