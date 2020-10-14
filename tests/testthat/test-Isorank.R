context("IsoRank")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
startm <- matrix(0, 10, 10)
diag(startm)[1:4] <- 1
seeds<-1:4

test_that("order of nodes getting matched", {
  expect_equal(graph_match_IsoRank(A, B, startm, seeds, method = "greedy")$order,
               c(1,2,3,4,6,10,8,7,5,9))
})
test_that("test LAP method", {
  expect_equal(graph_match_IsoRank(A, B, startm, seeds, method = "LAP")$ns, 4)
})

# create similarity score matrix that is not a square matrix
startm <- startm[1:8,]
test_that("test padding for similarity scores", {
  expect_equal(graph_match_IsoRank(A, B, startm, seeds = NULL, method = "LAP")$corr,
               data.frame(corr_A = c(1:10), corr_B = c(1,2,3,4,9,6,7,8,5,10), row.names = as.character(1:10)))
})

set.seed(12)
g2 <- sample_correlated_gnp_pair(n = 10, corr = 0.7, p = 0.3)
A2 <- g2$graph1
B2 <- g2$graph2

A_l <- list(A, A2)
B_l <- list(B, B2)

test_that("IsoRank multi-layer", {
  expect_equal(graph_match_IsoRank(A_l, B_l, startm, seeds = 1:3)$ns, 3)
})
