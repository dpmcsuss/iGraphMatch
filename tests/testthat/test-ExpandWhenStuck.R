context("ExpandWhenStuck")

set.seed(1)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
A <- cgnp_pair$graph1
B <- cgnp_pair$graph2

seeds <- 1:5
test_that("matching correspondence between graph1 and graph2", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds, r = 2)$corr,
               data.frame(corr_A = c(1:10), corr_B = c(1:5,8,10,7,9,6),row.names=as.character(1:10)))
})
test_that("order of nodes getting matched", {
  expect_equal(graph_match_ExpandWhenStuck(A, B, seeds, r = 2)$order,
               c(1,2,3,4,5,7,8,9,6,10))
})
test_that("test number of seeds", {
  expect_equal(graph_match_ExpandWhenStuck(A, B,seeds, r = 2)$ns, 5)
})
