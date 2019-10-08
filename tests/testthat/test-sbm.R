context("sbm & sbm with junk vertices")


#Sample graphs pair from stochastic block model
set.seed(123)
pm <- cbind( c(.1, .001), c(.001, .05) )
g<-sample_correlated_sbm_pair(10, pref.matrix=pm, block.sizes=c(3,7), rho=0.5,permutation=1:10,directed=TRUE,loops=TRUE)

test_that("number of vertices", {
  expect_equal(igraph::vcount(g$graph1),10) 
  expect_equal(igraph::vcount(g$graph2),10)
})
test_that("number of edges", {
  expect_equal(igraph::ecount(g$graph1),2)
  expect_equal(igraph::ecount(g$graph2),2)
})
test_that("degree of vertex in each graph", {
  expect_equal(igraph::degree(g$graph1),c(0,0,0,1,1,2,0,0,0,0))
  expect_equal(igraph::degree(g$graph2),c(0,0,0,1,1,1,0,0,0,1)) 
})


set.seed(123)
#Sample graphs pair from stochastic block model with junk vertices
pm <- cbind( c(.1, .001), c(.001, .05) )
g1<-sample_correlated_sbm_pair_w_junk(10, pref.matrix=pm, block.sizes=c(3,7), rho=0.5,core.block.sizes=c(2,5),permutation=1:10,directed=TRUE,loops=TRUE)

test_that("number of vertices", {
  expect_equal(igraph::vcount(g1$graph1),10) 
  expect_equal(igraph::vcount(g1$graph2),10)
})
test_that("number of edges", {
  expect_equal(igraph::ecount(g1$graph1),2)
  expect_equal(igraph::ecount(g1$graph2),1)
})
test_that("degree of vertex in each graph", {
  expect_equal(igraph::degree(g1$graph1),c(0,0,1,0,1,0,0,0,0,2))
  expect_equal(igraph::degree(g1$graph2),c(0,0,1,0,0,0,0,0,0,1)) 
})



