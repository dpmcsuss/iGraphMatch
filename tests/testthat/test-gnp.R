context("gnp & gnp with junk vertices")


#sample correlated G(n,p) random graphs
set.seed(123)
g<-sample_correlated_gnp_pair(10, 1, 0.5,directed=TRUE,loops=TRUE,permutation=c(sample(1:10)))

test_that("number of vertices", {
  expect_equal(igraph::vcount(g$graph1),10) 
  expect_equal(igraph::vcount(g$graph2),10)
})
test_that("number of edges", {
  expect_equal(igraph::ecount(g$graph1),52)
  expect_equal(igraph::ecount(g$graph2),52)
})
test_that("degree of vertex in each graph", {
  expect_equal(igraph::degree(g$graph1),c(11,9,7,9,10,13,11,9,13,12))
  expect_equal(igraph::degree(g$graph2),c(9,9,13,11,12,7,10,9,13,11)) 
})


set.seed(123)
#sample correlated G(n,p) random graphs with junk vertices
g<-sample_correlated_gnp_pair_w_junk(10, 1, 0.5, 4,directed=TRUE, loops=TRUE)

test_that("number of vertices", {
  expect_equal(igraph::vcount(g$graph1),10) 
  expect_equal(igraph::vcount(g$graph2),10)
})
test_that("number of edges", {
  expect_equal(igraph::ecount(g$graph1),51)
  expect_equal(igraph::ecount(g$graph2),54)
})
test_that("degree of vertex in each graph", {
  expect_equal(igraph::degree(g$graph1),c(10,12,10,7,11,11,8,14,11,8))
  expect_equal(igraph::degree(g$graph2),c(5,12,8,12,14,14,13,10,9,11)) 
})



