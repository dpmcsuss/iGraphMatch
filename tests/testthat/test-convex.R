context("Convex")

# sample pair of graphs w. 10 vertices
set.seed(123)
cgnp_pair <- sample_correlated_gnp_pair(n=10, corr=0.3, p=0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2
seeds <- 1:10 <= 3

test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_convex(A, B, seeds = seeds)$corr,
               data.frame(corr_A = c(1:10), corr_B = c(1,2,3,6,5,8,10,9,7,4)))
})
test_that("permutation matrix", {
  expect_equal(graph_match_convex(A, B)$P,get_perm(10, 10, data.frame(corr_A = c(1:10), corr_B = c(1,2,3,6,5,8,10,9,7,4))))
})
test_that("number of seeds", {
  expect_equal(graph_match_convex(A, B)$ns,3)
})

